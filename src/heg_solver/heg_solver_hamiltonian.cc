#include "heg_solver.h"

#include "../array_math.h"
#include "../det/det.h"
#include "../det/spin_det.h"

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
    for (const int p : occ_pq_up) H += sum(square(k_points[p] * k_unit) * 0.5);
    for (const int p : occ_pq_dn) H += sum(square(k_points[p] * k_unit) * 0.5);

    // Two electrons operator.
    for (std::size_t i = 0; i < n_up; i++) {
      const int p = occ_pq_up[i];
      for (std::size_t j = i + 1; j < n_up; j++) {
        const int q = occ_pq_up[j];
        H -= H_unit / sum(square(k_points[p] - k_points[q]));
      }
    }
    for (std::size_t i = 0; i < n_dn; i++) {
      const int p = occ_pq_dn[i];
      for (std::size_t j = i + 1; j < n_dn; j++) {
        const int q = occ_pq_dn[j];
        H -= H_unit / sum(square(k_points[p] - k_points[q]));
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

    H = H_unit / sum(square(k_points[orb_p] - k_points[orb_r]));
    if (n_eor_up != 2) H -= H_unit / sum(square(k_points[orb_p] - k_points[orb_s]));

    const int gamma_exp =
        get_gamma_exp(det_pq.up, eor_up_set_bits) + get_gamma_exp(det_pq.dn, eor_dn_set_bits) +
        get_gamma_exp(det_rs.up, eor_up_set_bits) + get_gamma_exp(det_rs.dn, eor_dn_set_bits);
    if ((gamma_exp & 1) == 1) H = -H;
  }
  return H;
}