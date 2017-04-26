#include "solver.h"

#include <boost/functional/hash.hpp>

#include "../det/det.h"
#include "../parallel.h"

void Solver::variation() {
  // Setup HF or existing wf as initial wf and evaluate energy.
  // Loop:
  //   Find connected determinants and add to wf.
  //   Evaluate new energy.
  //   Exit loop if energy change within threshold.

  const double THRESHOLD = 1.0e-6;

  if (wf.size() == 0) {
    Det det;
    for (int i = 0; i < n_up; i++) det.up.set_orb(i, true);
    for (int i = 0; i < n_dn; i++) det.dn.set_orb(i, true);
    wf.append_term(det, 1.0);
    energy_hf = energy_var = hamiltonian(det, det);
    printf("HF energy: %.10f\n", energy_hf);
  }

  double energy_var_new = 0.0;
  std::unordered_set<
      std::pair<Orbitals, Orbitals>,
      boost::hash<std::pair<Orbitals, Orbitals>>>
      var_dets_set;
  for (const auto& term : wf.get_terms()) var_dets_set.insert(term.first.encode());
  int cnt = 0;
  while (fabs(energy_var - energy_var_new) > THRESHOLD) {
    if (Parallel::get_id() == 0) printf("Variation Iteration: %d\n", ++cnt);
    energy_var = energy_var_new;
    std::list<Det> new_dets;
    for (const auto& term : wf.get_terms()) {
      const auto& connected_dets =
          find_connected_dets(term.first, fabs(eps_var / term.second));
      for (const auto& new_det : connected_dets) {
        if (var_dets_set.count(new_det.encode()) == 0) {
          var_dets_set.insert(new_det.encode());
          new_dets.push_back(new_det);
        }
      }
    }
    for (const auto& det : new_dets) wf.append_term(det, 0.0);
    new_dets.clear();
    if (Parallel::get_id() == 0) {
      printf("Number of dets: %d\n", static_cast<int>(var_dets_set.size()));
    }
  }
}