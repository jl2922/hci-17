#include "heg_solver.h"

#include <boost/format.hpp>
#include <boost/range/adaptor/reversed.hpp>

#include "../array_math.h"
#include "../big_unordered_map.h"
#include "../config.h"
#include "../parallel.h"
#include "../regression/linear_regression.h"
#include "../time/time.h"
#include "k_points_util.h"

void HEGSolver::solve() {
  printf("Proc %d running on %s\n", Parallel::get_id(), Parallel::get_host().c_str());
  n_up = Config::get<std::size_t>("n_up");
  n_dn = Config::get<std::size_t>("n_dn");
  rcut_vars = Config::get_array<double>("rcut_vars");
  eps_vars = Config::get_array<double>("eps_vars");
  rcut_pts = Config::get_array<double>("rcut_pts");
  eps_pts = Config::get_array<double>("eps_pts");
  std::sort(rcut_vars.begin(), rcut_vars.end());
  std::sort(eps_vars.begin(), eps_vars.end());
  std::sort(rcut_pts.begin(), rcut_pts.end());
  std::sort(eps_pts.begin(), eps_pts.end());
  std::reverse(eps_vars.begin(), eps_vars.end());
  std::reverse(eps_pts.begin(), eps_pts.end());

  Time::start("variation stage");
  for (const double rcut_var : rcut_vars) {
    std::string rcut_var_event = str(boost::format("variation with rcut_var: %#.4g") % rcut_var);
    Time::start(rcut_var_event);
    this->rcut_var = rcut_var;
    setup();
    for (const double eps_var : eps_vars) {
      std::string eps_var_event = str(boost::format("variation with eps_var: %#.4g") % eps_var);
      Time::start(eps_var_event);
      this->eps_var = eps_var;
      if (!load_variation_result()) {
        variation();
        save_variation_result();
      }
      Time::end(eps_var_event);
    }
    Time::end(rcut_var_event);
  }
  Time::end("variation stage");

  Time::start("perturbation stage");
  n_orbs_pts.clear();
  for (const double rcut_pt : rcut_pts) {
    n_orbs_pts.push_back(KPointsUtil::get_n_k_points(rcut_pt) * 2);
  }
  // Start from the largest PT so that it fails earlier upon insufficient memory.
  for (const double rcut_var : rcut_vars | boost::adaptors::reversed) {
    std::string rcut_var_event = str(boost::format("perturbation with rcut_var: %#.4g") % rcut_var);
    Time::start(rcut_var_event);
    this->rcut_var = rcut_var;
    for (const double eps_var : eps_vars | boost::adaptors::reversed) {
      std::string eps_var_event = str(boost::format("perturbation with eps_var: %#.4g") % eps_var);
      Time::start(eps_var_event);
      this->eps_var = eps_var;
      assert(load_variation_result());
      perturbation();
      Time::end(eps_var_event);
    }
    Time::end(rcut_var_event);
  }
  Time::end("perturbation stage");

  Time::start("extrapolation");
  extrapolate();
  Time::end("extrapolation");
}

void HEGSolver::setup() {
  const double r_s = Config::get<double>("r_s");
  const double density = 3.0 / (4.0 * M_PI * pow(r_s, 3));
  const double cell_length = pow((n_up + n_dn) / density, 1.0 / 3);
  k_unit = 2 * M_PI / cell_length;
  H_unit = 1.0 / (M_PI * cell_length);

  k_points = KPointsUtil::generate_k_points(rcut_var);
  k_lut = KPointsUtil::generate_k_lut(k_points);

  Time::start("Generate HCI queue.");
  generate_hci_queue(rcut_var);
  Time::end("Generate HCI queue.");

  wf.clear();
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
      if (sum(square(diff_pr)) == sum(square(diff_ps))) continue;
      const double abs_H = fabs(1.0 / sum(square(diff_pr)) - 1.0 / sum(square(diff_ps)));
      if (abs_H < DBL_EPSILON) continue;
      const auto& item = TinyInt3Double(cast<TinyInt>(diff_pr), abs_H * H_unit);
      same_spin_hci_queue[cast<TinyInt>(diff_pq)].push_back(item);
    }
  }
  for (auto& kv : same_spin_hci_queue) {
    auto& items = kv.second;
    std::stable_sort(
        items.begin(), items.end(), [](const TinyInt3Double& a, const TinyInt3Double& b) -> bool {
          return a.second > b.second;
        });
    max_abs_H = std::max(max_abs_H, items.front().second);
  }

  // Opposite spin.
  for (const auto& diff_pr : k_diffs) {
    const double abs_H = 1.0 / sum(square(diff_pr));
    if (abs_H < DBL_EPSILON) continue;
    const auto& item = TinyInt3Double(cast<TinyInt>(diff_pr), abs_H * H_unit);
    opposite_spin_hci_queue.push_back(item);
  }
  std::stable_sort(
      opposite_spin_hci_queue.begin(),
      opposite_spin_hci_queue.end(),
      [](const TinyInt3Double& a, const TinyInt3Double& b) -> bool { return a.second > b.second; });
  max_abs_H = std::max(max_abs_H, opposite_spin_hci_queue.front().second);
}

void HEGSolver::save_variation_result() {
  if (Parallel::get_id() != 0) return;
  std::ofstream var_file;
  std::string filename = str(boost::format("var_%.5f_%.3f.txt") % eps_var % rcut_var);
  var_file.open(filename);
  var_file << boost::format("%.15g %.15g\n") % energy_hf % energy_var;
  var_file << boost::format("%d %d %d\n") % n_up % n_dn % wf.size();
  for (const auto& term : wf.get_terms()) {
    var_file << boost::format("%.15g\n") % term.coef;
    var_file << term.det.up << std::endl << term.det.dn << std::endl;
  }
  var_file.close();
  printf("Variation result saved to: %s\n", filename.c_str());
}

bool HEGSolver::load_variation_result() {
  std::ifstream var_file;
  std::string filename = str(boost::format("var_%.5f_%.3f.txt") % eps_var % rcut_var);
  std::size_t wf_size;
  int orb_id;
  double coef;
  var_file.open(filename);
  if (!var_file.is_open()) return false;  // Does not exist.
  var_file >> energy_hf >> energy_var;
  var_file >> n_up >> n_dn >> wf_size;
  wf.clear();
  for (std::size_t i = 0; i < wf_size; i++) {
    var_file >> coef;
    Det det;
    for (std::size_t j = 0; j < n_up; j++) {
      var_file >> orb_id;
      det.up.set_orb(orb_id, true);
    }
    for (std::size_t j = 0; j < n_dn; j++) {
      var_file >> orb_id;
      det.dn.set_orb(orb_id, true);
    }
    wf.append_term(det, coef);
  }
  var_file.close();
  if (Parallel::get_id() == 0)
    printf("Loaded %'d dets from: %s\n", static_cast<int>(wf_size), filename.c_str());
  return true;
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
    Int3 k_change;
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

std::list<IntPair> get_pq_pairs(const Det& det, const int dn_offset) {
  const auto& occ_up = det.up.get_elec_orbs();
  const auto& occ_dn = det.dn.get_elec_orbs();
  const std::size_t n_up = det.up.get_n_elecs();
  const std::size_t n_dn = det.dn.get_n_elecs();

  std::list<IntPair> pq_pairs;

  for (std::size_t i = 0; i < n_up; i++) {
    for (std::size_t j = i + 1; j < n_up; j++) {
      pq_pairs.push_back(IntPair(occ_up[i], occ_up[j]));
    }
  }
  for (std::size_t i = 0; i < n_dn; i++) {
    for (std::size_t j = i + 1; j < n_dn; j++) {
      pq_pairs.push_back(IntPair(occ_dn[i] + dn_offset, occ_dn[j] + dn_offset));
    }
  }
  for (std::size_t i = 0; i < n_up; i++) {
    for (std::size_t j = 0; j < n_dn; j++) {
      pq_pairs.push_back(IntPair(occ_up[i], occ_dn[j] + dn_offset));
    }
  }

  return pq_pairs;
}

std::list<Det> HEGSolver::find_connected_dets(const Det& det, const double eps) const {
  std::list<Det> connected_dets;
  connected_dets.push_back(det);

  if (max_abs_H < eps) return connected_dets;

  const int dn_offset = static_cast<int>(k_points.size());
  const auto& pq_pairs = get_pq_pairs(det, dn_offset);

  for (const auto& pq_pair : pq_pairs) {
    const int p = pq_pair.first;
    const int q = pq_pair.second;

    // Get rs pairs.
    int pp = p, qq = q;
    if (p >= dn_offset && q >= dn_offset) {
      pp -= dn_offset;
      qq -= dn_offset;
    } else if (p < dn_offset && q >= dn_offset && p > q - dn_offset) {
      pp = q - dn_offset;
      qq = p + dn_offset;
    }
    bool same_spin = false;
    std::vector<TinyInt3Double> const* items_ptr;
    if (pp < dn_offset && qq < dn_offset) {
      same_spin = true;
      const auto& diff_pq = cast<TinyInt>(k_points[qq] - k_points[pp]);
      items_ptr = &(same_spin_hci_queue.find(diff_pq)->second);
    } else {
      items_ptr = &(opposite_spin_hci_queue);
    }
    const auto& items = *items_ptr;
    int qs_offset = 0;
    if (!same_spin) qs_offset = dn_offset;

    for (const auto& item : items) {
      if (item.second < eps) break;
      const auto& diff_pr = cast<int>(item.first);
      const auto it_r = k_lut.find(diff_pr + k_points[pp]);
      if (it_r == k_lut.end()) continue;
      int r = it_r->second;
      const auto it_s = k_lut.find(k_points[pp] + k_points[qq - qs_offset] - k_points[r]);
      if (it_s == k_lut.end()) continue;
      int s = it_s->second;
      if (same_spin && s < r) continue;
      s += qs_offset;
      if (p >= dn_offset && q >= dn_offset) {
        r += dn_offset;
        s += dn_offset;
      } else if (p < dn_offset && q >= dn_offset && p > q - dn_offset) {
        const int tmp = s;
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

#ifndef SERIAL
template <>
std::size_t BigUnorderedMap<PTKey, double, boost::hash<PTKey>>::get_target(const PTKey& key) {
  boost::hash<OrbitalsPair> det_code_hasher;
  return proc_map[det_code_hasher(key.first) % total_proc_buckets];
}
#endif

void HEGSolver::perturbation() {
  // Perform perturbation with smallest eps and largest rcut.
  const double rcut_pt_max = rcut_pts.back();
  const double eps_pt_min = eps_pts.back();
  k_points = KPointsUtil::generate_k_points(rcut_pt_max);
  k_lut = KPointsUtil::generate_k_lut(k_points);
  generate_hci_queue(rcut_pt_max);
  if (Parallel::get_id() == 0) {
    printf("PT with rcut_pt_max = %#.4g, eps_pt_min = %#.4g\n", rcut_pt_max, eps_pt_min);
    printf("Number of max perturbation orbitals: %d\n", static_cast<int>(k_points.size() * 2));
  }

  // Cache variation determinants.
  var_dets_set.clear();
  const auto& terms = wf.get_terms();
  var_dets_set.rehash(terms.size() * 2);  // <20% conflict rate with 50% hash load.
  for (const auto& term : terms) {
    var_dets_set.insert(term.det.encode());
  }

  Time::start("setup hash table");
  unsigned long long n_pt_dets_estimate = estimate_n_pt_dets(eps_pt_min);
  if (Parallel::get_id() == 0) printf("Estimated PT terms: %'llu\n", n_pt_dets_estimate);
  std::pair<PTKey, double> skeleton;  // For reducing the amount of MPI data transfer.
  skeleton.first.first = wf.get_terms().front().det.encode();
  BigUnorderedMap<PTKey, double, boost::hash<PTKey>> pt_sums(skeleton);
  pt_sums.reserve(n_pt_dets_estimate * 2);
  unsigned long long hash_buckets = pt_sums.bucket_count();
  if (Parallel::get_id() == 0) printf("Reserved %'llu total hash buckets.\n", hash_buckets);
  Time::end("setup hash table");

  Time::start("search for perturbation dets");
  int progress = 1;  // For print.
  std::size_t i = 0;
  const std::size_t n = wf.size();
  for (const auto& term : wf.get_terms()) {
    if ((i++) % Parallel::get_n() != static_cast<std::size_t>(Parallel::get_id())) continue;
    const auto& connected_dets = find_connected_dets(term.det, eps_pt_min / fabs(term.coef));
    for (const auto& det_a : connected_dets) {
      if (var_dets_set.count(det_a.encode()) == 1) continue;
      const double H_ai = hamiltonian(term.det, det_a);
      if (fabs(H_ai) < DBL_EPSILON) continue;
      const double partial_sum = H_ai * term.coef;
      PTCategory category = get_pt_category((det_a.get_highest_orb() + 1) * 2, fabs(partial_sum));
      PTKey ptKey(det_a.encode(), category);
      pt_sums.async_inc(ptKey, partial_sum);
    }
    if (Parallel::get_id() == 0 && i >= n / 100 * progress) {
      const auto& local_map = pt_sums.get_local_map();
      Time::checkpoint("search for perturbation dets");
      printf(
          "Master progress: %d%%. Local PT keys: %'lu, hash load: %.2f\n",
          progress,
          local_map.size(),
          local_map.load_factor());
      progress *= 2;
    }
  }
  pt_sums.complete_async_incs();
  unsigned long long n_pt_keys = pt_sums.size();
  if (Parallel::get_id() == 0) printf("Total PT keys: %'llu\n", n_pt_keys);
  Time::end("search for perturbation dets");

  Time::start("accumulate contributions");
  const auto& local_map = pt_sums.get_local_map();
  for (const double rcut_pt : rcut_pts) {
    std::size_t n_orbs_pt = KPointsUtil::get_n_k_points(rcut_pt) * 2;
    std::string n_orbs_pt_event = "accumulate for n_orbs_pt: " + std::to_string(n_orbs_pt);
    Time::start(n_orbs_pt_event);
    for (const double eps_pt : eps_pts) {
      std::string eps_pt_event = str(boost::format("accumulate for eps_pt: %.4g") % eps_pt);
      Time::start(eps_pt_event);
      const auto& related_categories = get_related_pt_categories(n_orbs_pt, eps_pt);
      energy_pt = 0.0;
      BigUnsignedInt n_pt_dets = 0;
      for (const auto& kv : local_map) {
        const auto& key = kv.first;
        const PTCategory category = key.second;
        if (std::find(related_categories.begin(), related_categories.end(), category) ==
            related_categories.end()) {
          continue;
        }
        bool is_smallest = true;
        double partial_sum = 0.0;
        for (const auto related_category : related_categories) {
          const PTKey related_key(key.first, related_category);
          if (local_map.count(related_key) == 1) {
            // Only the smallest one submits the contribution.
            if (related_category < category) {
              is_smallest = false;
              break;
            }
            partial_sum += local_map.at(related_key);
          }
        }
        if (is_smallest) {
          Det det_a;
          det_a.decode(key.first);
          const double H_aa = hamiltonian(det_a, det_a);
          energy_pt += pow(partial_sum, 2) / (energy_var - H_aa);
          n_pt_dets++;
        }
      }  // local_map loop.
      Parallel::reduce_to_sum(n_pt_dets);  // DEBUG.
      Parallel::reduce_to_sum(energy_pt);
      const double correlation_energy = energy_var + energy_pt - energy_hf;
      std::size_t n_orbs_var = KPointsUtil::get_n_k_points(rcut_var) * 2;
      if (Parallel::get_id() == 0) {
        printf("Number of related PT dets: %'llu\n", n_pt_dets);
        printf("n_orbs_var: %d\n", static_cast<int>(n_orbs_var));
        printf("eps_var: %#.4g\n", eps_var);
        printf("n_orbs_pt: %d\n", static_cast<int>(n_orbs_pt));
        printf("eps_pt: %#.4g\n", eps_pt);
        printf("Perturbation energy: %#.12g Ha\n", energy_pt);
        printf("Correlation Energy: %.12g Ha\n", correlation_energy);
        std::vector<double> parameter_set({1.0 / n_orbs_var, eps_var, 1.0 / n_orbs_pt, eps_pt});
        parameter_sets.push_back(parameter_set);
        printf("Number of parameter sets: %d\n", static_cast<int>(parameter_sets.size()));
        results.push_back(correlation_energy);
      }
      Time::end(eps_pt_event);
    }  // eps_pts loop.
    Time::end(n_orbs_pt_event);
  }  // n_orbs_pts loop.
  Time::end("accumulate contributions");
}

std::vector<PTCategory> HEGSolver::get_related_pt_categories(
    const std::size_t n_orbs, const double eps) {
  std::vector<PTCategory> related_categories;
  for (const std::size_t n_orbs_pt : n_orbs_pts) {
    if (n_orbs < n_orbs_pt) continue;
    for (const double eps_pt : eps_pts) {
      if (eps > eps_pt) continue;
      related_categories.push_back(get_pt_category(n_orbs_pt, eps_pt));
    }
  }
  return related_categories;
}

PTCategory HEGSolver::get_pt_category(const std::size_t n_orbs, const double eps) {
  // The higher 4 bits represent eps_pt categories and the lower 4 bits represent n_orbs_pt.
  PTCategory category_eps = 0;
  PTCategory category_n_orbs = 0;
  for (const double eps_pt : eps_pts) {
    if (eps >= eps_pt) category_eps++;
  }
  category_eps--;
  for (const std::size_t n_orbs_pt : n_orbs_pts) {
    if (n_orbs <= n_orbs_pt) category_n_orbs++;
  }
  category_n_orbs--;
  return category_eps + (category_n_orbs << 4);
}

void HEGSolver::extrapolate() {
  if (Parallel::get_id() == 0) {
    // Second order surface fit.
    std::vector<std::string> basic_quantities({"1/n_orbs_var", "eps_var", "1/n_orbs_pt", "eps_pt"});
    for (int i = 0; i < 4; i++) parameter_names.push_back(basic_quantities[i]);
    for (int i = 0; i < 4; i++) {
      for (int j = i; j < 4; j++) {
        for (auto& parameter_set : parameter_sets) {
          parameter_set.push_back(parameter_set[i] * parameter_set[j]);
        }
        parameter_names.push_back(parameter_names[i] + " * " + parameter_names[j]);
      }
    }
    parameter_names.push_back("Intercept");

    try {
      LinearRegression lr(parameter_sets, results);
      const auto& estimate = lr.get_estimate();
      const auto& stdev = lr.get_stdev();
      const auto& prob_t = lr.get_prob_t();
      printf("%30s %20s %15s %15s\n", "parameter", "estimate", "stdev", "P>|t|");
      for (int i = 0; i < 15; i++) {
        printf(
            "%30s %#20.10g %#15.5g %#15.5g\n",
            parameter_names[i].c_str(),
            estimate[i],
            stdev[i],
            prob_t[i]);
      }
    } catch (std::exception& e) {
      std::cout << e.what() << std::endl;
      return;
    }
  }
}
