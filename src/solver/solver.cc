#include "solver.h"

#include "../std.h"

Det Solver::generate_hf_det() {
  Det det;
  for (std::size_t i = 0; i < n_up; i++) det.up.set_orb(i, true);
  for (std::size_t i = 0; i < n_dn; i++) det.dn.set_orb(i, true);
  return det;
}

void Solver::variation(const double eps_var) {
  const double THRESHOLD = 1.0e-6;

  // Setup HF or existing wf as initial wf and evaluate energy.
  if (wf.size() == 0) {
    const Det& det_hf = generate_hf_det();
    wf.append_term(det_hf, 1.0);
    // energy_hf = energy_var = hamiltonian(det_hf, det_hf);
    // if (Parallel::get_id() == 0) printf("HF energy: %#.15g Ha\n", energy_hf);
  }

  // double energy_var_new = 0.0;  // Ensures the first iteration will run.
  // var_dets_set.clear();
  // for (const auto& term : wf.get_terms()) var_dets_set.insert(term.det.encode());
  // int iteration = 0;  // For print.
  // while (fabs(energy_var - energy_var_new) > THRESHOLD) {
  //   Time::start("Variation Iteration: " + std::to_string(iteration));
  //   energy_var = energy_var_new;

  //   // Find connected determinants.
  //   std::list<Det> new_dets;
  //   for (const auto& term : wf.get_terms()) {
  //     const auto& connected_dets = find_connected_dets(term.det, fabs(eps_var / term.coef));
  //     for (const auto& new_det : connected_dets) {
  //       if (var_dets_set.count(new_det.encode()) == 0) {
  //         var_dets_set.insert(new_det.encode());
  //         new_dets.push_back(new_det);
  //       }
  //     }
  //   }
  //   for (const auto& det : new_dets) wf.append_term(det, 0.0);
  //   new_dets.clear();
  //   if (Parallel::get_id() == 0) {
  //     printf("Number of dets: %'llu\n", static_cast<BigUnsignedInt>(var_dets_set.size()));
  //   }
  //   energy_var_new = diagonalize();
  //   if (Parallel::get_id() == 0) printf("Variation energy: %#.15g Ha\n", energy_var_new);

  //   Time::end("Variation Iteration: " + std::to_string(iteration));
  //   iteration++;
  // }

  // energy_var = energy_var_new;
  // if (Parallel::get_id() == 0) printf("Final variation energy: %#.15g Ha\n", energy_var);
}