#include "heg_solver.h"

#include "../array_math.h"
#include "../types.h"

std::list<IntPair> get_pq_pairs(const Det& det, const int dn_offset) {
  const auto& occ_up = det.up.get_elec_orbs();
  const auto& occ_dn = det.dn.get_elec_orbs();
  const int n_up = det.up.get_n_elecs();
  const int n_dn = det.dn.get_n_elecs();

  std::list<IntPair> pq_pairs;

  for (int i = 0; i < n_up; i++) {
    for (int j = i + 1; j < n_up; j++) {
      pq_pairs.push_back(IntPair(occ_up[i], occ_up[j]));
    }
  }
  for (int i = 0; i < n_dn; i++) {
    for (int j = i + 1; j < n_dn; j++) {
      pq_pairs.push_back(IntPair(occ_dn[i] + dn_offset, occ_dn[j] + dn_offset));
    }
  }
  for (int i = 0; i < n_up; i++) {
    for (int j = 0; j < n_dn; j++) {
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