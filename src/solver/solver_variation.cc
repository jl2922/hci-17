#include "solver.h"

#include <boost/functional/hash.hpp>

#include "../det/det.h"
#include "../parallel.h"
#include "../time.h"
#include "diagonalization/davidson.h"

Det generate_hf_det(const std::size_t n_up, const std::size_t n_dn) {
  Det det;
  for (std::size_t i = 0; i < n_up; i++) det.up.set_orb(i, true);
  for (std::size_t i = 0; i < n_dn; i++) det.dn.set_orb(i, true);
  return det;
}

void Solver::variation() {
  // Setup HF or existing wf as initial wf and evaluate energy.
  // Loop:
  //   Find connected determinants and add to wf.
  //   Evaluate new energy.
  //   Exit loop if energy change within threshold.

  const double THRESHOLD = 1.0e-6;

  if (wf.size() == 0) {
    const Det& det_hf = generate_hf_det(n_up, n_dn);
    wf.append_term(det_hf, 1.0);
    energy_hf = energy_var = hamiltonian(det_hf, det_hf);
    if (Parallel::get_id() == 0) printf("HF energy: %#.15g Ha\n", energy_hf);
  }

  double energy_var_new = 0.0;  // Ensures the first iteration will run.
  std::unordered_set<OrbitalsPair, boost::hash<OrbitalsPair>> var_dets_set;
  for (const auto& term : wf.get_terms()) var_dets_set.insert(term.det.encode());
  int iteration = 0;
  while (fabs(energy_var - energy_var_new) > THRESHOLD) {
    Time::start("Variation Iteration: " + std::to_string(iteration));
    energy_var = energy_var_new;

    // Find connected determinants.
    std::list<Det> new_dets;
    for (const auto& term : wf.get_terms()) {
      const auto& connected_dets = find_connected_dets(term.det, fabs(eps_var / term.coef));
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
      printf("Number of dets: %'llu\n", static_cast<BigUnsignedInt>(var_dets_set.size()));
    }
    energy_var_new = diagonalize();
    if (Parallel::get_id() == 0) printf("Variation energy: %#.15g Ha\n", energy_var_new);

    Time::end("Variation Iteration: " + std::to_string(iteration));
    iteration++;
  }

  energy_var = energy_var_new;
  if (Parallel::get_id() == 0) printf("Final variation energy: %#.15g Ha\n", energy_var);
}

double Solver::diagonalize() {
  std::vector<double> diagonal;
  std::vector<double> initial_vector;
  diagonal.reserve(wf.size());
  initial_vector.reserve(wf.size());
  for (const auto& term : wf.get_terms()) {
    const auto& det = term.det;
    diagonal.push_back(hamiltonian(det, det));
    initial_vector.push_back(term.coef);
  }
  HelperStrings helper_strings(wf.get_dets());
  std::function<std::vector<double>(std::vector<double>)> apply_hamiltonian_func =
      std::bind(&Solver::apply_hamiltonian, this, std::placeholders::_1, helper_strings);

  Davidson davidson(diagonal, apply_hamiltonian_func, wf.size());
  if (Parallel::get_id() == 0) davidson.set_verbose(true);
  davidson.diagonalize(initial_vector);

  double energy_var = davidson.get_lowest_eigenvalue();
  const auto& coefs_new = davidson.get_lowest_eigenvector();
  wf.set_coefs(coefs_new);
  wf.sort_by_coefs();

  return energy_var;
}