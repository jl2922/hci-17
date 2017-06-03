#include "heg_solver.h"

#include "../array_math.h"
#include "../config.h"
#include "../std.h"
#include "../timer.h"
#include "k_points_util.h"

void HEGSolver::solve() {
  // Get configurations.
  n_up = Config::get<std::size_t>("n_up");
  n_dn = Config::get<std::size_t>("n_dn");
  rcut_vars = Config::get_array<double>("rcut_vars");
  eps_vars = Config::get_array<double>("eps_vars");
  rcut_pts = Config::get_array<double>("rcut_pts");
  eps_pts = Config::get_array<double>("eps_pts");
  std::sort(rcut_vars.begin(), rcut_vars.end());
  std::sort(rcut_pts.begin(), rcut_pts.end());
  std::sort(eps_vars.begin(), eps_vars.end());
  std::sort(eps_pts.begin(), eps_pts.end());
  std::reverse(eps_vars.begin(), eps_vars.end());
  std::reverse(eps_pts.begin(), eps_pts.end());

  // Variation.
  Timer::start("variation");
  for (const double rcut_var : rcut_vars) {
    Timer::start(str(boost::format("rcut_var: %#.4g") % rcut_var));
    Timer::start("setup");
    setup(rcut_var);
    Timer::end();
    for (const double eps_var : eps_vars) {
      Timer::start(str(boost::format("eps_var: %#.4g") % eps_var));
      variation(eps_var);
      Timer::end();
    }
    Timer::end();
  }
  Timer::end();

  // Perturbation.

  // Extrapolation.
}

void HEGSolver::setup(const double rcut) {
  const double r_s = Config::get<double>("r_s");
  const double density = 3.0 / (4.0 * M_PI * pow(r_s, 3));
  const double cell_length = pow((n_up + n_dn) / density, 1.0 / 3);
  k_unit = 2 * M_PI / cell_length;
  H_unit = 1.0 / (M_PI * cell_length);

  k_points = KPointsUtil::generate_k_points(rcut);
  k_lut = KPointsUtil::generate_k_lut(k_points);

  generate_hci_queue(rcut);
}

void HEGSolver::generate_hci_queue(const double rcut) {
  same_spin_hci_queue.clear();
  opposite_spin_hci_queue.clear();
  max_abs_H = 0.0;

  // Common dependencies.
  const auto& k_diffs = KPointsUtil::get_k_diffs(k_points);

  // Same spin.
  for (const auto& diff_pq : k_diffs) {
    for (const auto& diff_pr : k_diffs) {
      const auto& diff_sr = diff_pr + diff_pr - diff_pq;  // Momentum conservation.
      if (diff_sr == 0 || norm(diff_sr) > rcut * 2) continue;
      const auto& diff_ps = diff_pr - diff_sr;
      if (diff_ps == 0) continue;
      if (squared_norm(diff_pr) == squared_norm(diff_ps)) continue;
      const double abs_H = fabs(1.0 / squared_norm(diff_pr) - 1.0 / squared_norm(diff_ps));
      if (abs_H < DBL_EPSILON) continue;
      const auto& item = KPointDouble(diff_pr, abs_H * H_unit);
      same_spin_hci_queue[diff_pq].push_back(item);
    }
  }
  for (auto& kv : same_spin_hci_queue) {
    auto& items = kv.second;
    std::stable_sort(
        items.begin(), items.end(), [](const KPointDouble& a, const KPointDouble& b) -> bool {
          return a.second > b.second;
        });
    max_abs_H = std::max(max_abs_H, items.front().second);
  }

  // Opposite spin.
  for (const auto& diff_pr : k_diffs) {
    const double abs_H = 1.0 / squared_norm(diff_pr);
    if (abs_H < DBL_EPSILON) continue;
    const auto& item = KPointDouble(diff_pr, abs_H * H_unit);
    opposite_spin_hci_queue.push_back(item);
  }
  std::stable_sort(
      opposite_spin_hci_queue.begin(),
      opposite_spin_hci_queue.end(),
      [](const KPointDouble& a, const KPointDouble& b) -> bool { return a.second > b.second; });
  max_abs_H = std::max(max_abs_H, opposite_spin_hci_queue.front().second);
}

int get_gamma_exp(const SpinDet& spin_det, const std::vector<Orbital>& eor) {
  int gamma_exp = 0;
  int ptr = 0;
  const auto& occ = spin_det.get_elec_orbs();
  for (const Orbital orb_id : eor) {
    if (!spin_det.get_orb(orb_id)) continue;
    while (occ[ptr] < orb_id) ptr++;
    gamma_exp += ptr;
  }
  return gamma_exp;
}

double HEGSolver::hamiltonian(const Det& det_pq, const Det& det_rs) const {
  double H = 0.0;

  if (det_pq == det_rs) {
    const auto& occ_pq_up = det_pq.up.get_elec_orbs();
    const auto& occ_pq_dn = det_pq.dn.get_elec_orbs();

    // One electron operator.
    for (const int p : occ_pq_up) H += squared_norm(k_points[p] * k_unit) * 0.5;
    for (const int p : occ_pq_dn) H += squared_norm(k_points[p] * k_unit) * 0.5;

    // Two electrons operator.
    for (std::size_t i = 0; i < n_up; i++) {
      const int p = occ_pq_up[i];
      for (std::size_t j = i + 1; j < n_up; j++) {
        const int q = occ_pq_up[j];
        H -= H_unit / squared_norm(k_points[p] - k_points[q]);
      }
    }
    for (std::size_t i = 0; i < n_dn; i++) {
      const int p = occ_pq_dn[i];
      for (std::size_t j = i + 1; j < n_dn; j++) {
        const int q = occ_pq_dn[j];
        H -= H_unit / squared_norm(k_points[p] - k_points[q]);
      }
    }
  } else {
    // Off-diagonal elements.
    Det det_eor;
    det_eor.from_eor(det_pq, det_rs);
    const std::size_t n_eor_up = det_eor.up.get_n_elecs();
    const std::size_t n_eor_dn = det_eor.dn.get_n_elecs();
    if (n_eor_up + n_eor_dn != 4) return 0.0;
    const auto& eor_up_set_bits = det_eor.up.get_elec_orbs();
    const auto& eor_dn_set_bits = det_eor.dn.get_elec_orbs();
    bool k_p_set = false, k_r_set = false;
    int orb_p = 0, orb_r = 0, orb_s = 0;

    // Obtain p, q, s.
    KPoint k_change;
    k_change.fill(0);
    for (const int orb_i : eor_up_set_bits) {
      if (det_pq.up.get_orb(orb_i)) {
        k_change -= k_points[orb_i];
        if (!k_p_set) {
          orb_p = orb_i;
          k_p_set = true;
        }
      } else {
        k_change += k_points[orb_i];
        if (!k_r_set) {
          orb_r = orb_i;
          k_r_set = true;
        } else {
          orb_s = orb_i;
        }
      }
    }
    for (const int orb_i : eor_dn_set_bits) {
      if (det_pq.dn.get_orb(orb_i)) {
        k_change -= k_points[orb_i];
        if (!k_p_set) {
          orb_p = orb_i;
          k_p_set = true;
        }
      } else {
        k_change += k_points[orb_i];
        if (!k_r_set) {
          orb_r = orb_i;
          k_r_set = true;
        } else {
          orb_s = orb_i;
        }
      }
    }

    // Check for momentum conservation.
    if (k_change != 0) return 0.0;

    H = H_unit / squared_norm(k_points[orb_p] - k_points[orb_r]);
    if (n_eor_up != 2) H -= H_unit / squared_norm(k_points[orb_p] - k_points[orb_s]);

    const int gamma_exp =
        get_gamma_exp(det_pq.up, eor_up_set_bits) + get_gamma_exp(det_pq.dn, eor_dn_set_bits) +
        get_gamma_exp(det_rs.up, eor_up_set_bits) + get_gamma_exp(det_rs.dn, eor_dn_set_bits);
    if ((gamma_exp & 1) == 1) H = -H;
  }
  return H;
}

std::list<Det> HEGSolver::find_connected_dets(const Det& det, const double eps) const {
  std::list<Det> connected_dets;
  connected_dets.push_back(det);

  if (max_abs_H < eps) return connected_dets;

  const std::size_t dn_offset = k_points.size();
  const auto& pq_pairs = get_pq_pairs(det, dn_offset);

  for (const auto& pq_pair : pq_pairs) {
    const std::size_t p = pq_pair.first;
    const std::size_t q = pq_pair.second;

    // Get rs pairs.
    std::size_t pp = p, qq = q;
    if (p >= dn_offset && q >= dn_offset) {
      pp -= dn_offset;
      qq -= dn_offset;
    } else if (p < dn_offset && q >= dn_offset && p > q - dn_offset) {
      pp = q - dn_offset;
      qq = p + dn_offset;
    }
    bool same_spin = false;
    std::vector<KPointDouble> const* items_ptr;
    if (pp < dn_offset && qq < dn_offset) {
      same_spin = true;
      const auto& diff_pq = k_points[qq] - k_points[pp];
      items_ptr = &(same_spin_hci_queue.find(diff_pq)->second);
    } else {
      items_ptr = &(opposite_spin_hci_queue);
    }
    const auto& items = *items_ptr;
    std::size_t qs_offset = 0;
    if (!same_spin) qs_offset = dn_offset;

    for (const auto& item : items) {
      if (item.second < eps) break;
      const auto& diff_pr = item.first;
      const auto it_r = k_lut.find(diff_pr + k_points[pp]);
      if (it_r == k_lut.end()) continue;
      std::size_t r = it_r->second;
      const auto it_s = k_lut.find(k_points[pp] + k_points[qq - qs_offset] - k_points[r]);
      if (it_s == k_lut.end()) continue;
      std::size_t s = it_s->second;
      if (same_spin && s < r) continue;
      s += qs_offset;
      if (p >= dn_offset && q >= dn_offset) {
        r += dn_offset;
        s += dn_offset;
      } else if (p < dn_offset && q >= dn_offset && p > q - dn_offset) {
        const std::size_t tmp = s;
        s = r + dn_offset;
        r = tmp - dn_offset;
      }

      // Test whether pqrs is a valid excitation for det.
      if (det.get_orb(r, dn_offset) || det.get_orb(s, dn_offset)) continue;
      connected_dets.push_back(det);
      Det& new_det = connected_dets.back();
      new_det.set_orb(p, dn_offset, false);
      new_det.set_orb(q, dn_offset, false);
      new_det.set_orb(r, dn_offset, true);
      new_det.set_orb(s, dn_offset, true);
    }
  }

  return connected_dets;
};
