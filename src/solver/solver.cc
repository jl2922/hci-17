#include "solver.h"

#include <boost/functional/hash.hpp>
#include "../std.h"
#include "../timer.h"
#include "../wavefunction/wavefunction.h"

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
    energy_hf = energy_var = hamiltonian(det_hf, det_hf);
    // if (Parallel::get_id() == 0) printf("HF energy: %#.15g Ha\n", energy_hf);
    printf("HF energy: %#.15g Ha\n", energy_hf);
  }

  double energy_var_new = 0.0;  // Ensures the first iteration will run.
  std::unordered_set<OrbitalsPair, boost::hash<OrbitalsPair>> var_dets_set;
  var_dets_set.clear();
  for (const auto& term : wf.get_terms()) var_dets_set.insert(term.det.encode());
  int iteration = 0;  // For print.
  while (fabs(energy_var - energy_var_new) > THRESHOLD) {
    Timer::start("Iteration: " + std::to_string(iteration));
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
    printf("Number of dets: %'llu\n", static_cast<unsigned long long>(var_dets_set.size()));
    // if (Parallel::get_id() == 0) {
    //   printf("Number of dets: %'llu\n", static_cast<BigUnsignedInt>(var_dets_set.size()));
    // }
    // energy_var_new = diagonalize();
    // if (Parallel::get_id() == 0) printf("Variation energy: %#.15g Ha\n", energy_var_new);

    Timer::end();
    iteration++;
  }

  energy_var = energy_var_new;
  // if (Parallel::get_id() == 0) printf("Final variation energy: %#.15g Ha\n", energy_var);
}

std::list<PQPair> Solver::get_pq_pairs(const Det& det, const std::size_t dn_offset) const {
  const auto& occ_up = det.up.get_elec_orbs();
  const auto& occ_dn = det.dn.get_elec_orbs();
  const std::size_t n_up = det.up.get_n_elecs();
  const std::size_t n_dn = det.dn.get_n_elecs();

  std::list<PQPair> pq_pairs;

  for (std::size_t i = 0; i < n_up; i++) {
    for (std::size_t j = i + 1; j < n_up; j++) {
      pq_pairs.push_back(PQPair(occ_up[i], occ_up[j]));
    }
  }
  for (std::size_t i = 0; i < n_dn; i++) {
    for (std::size_t j = i + 1; j < n_dn; j++) {
      pq_pairs.push_back(PQPair(occ_dn[i] + dn_offset, occ_dn[j] + dn_offset));
    }
  }
  for (std::size_t i = 0; i < n_up; i++) {
    for (std::size_t j = 0; j < n_dn; j++) {
      pq_pairs.push_back(PQPair(occ_up[i], occ_dn[j] + dn_offset));
    }
  }

  return pq_pairs;
}

// std::vector<double> Solver::diagonalize() {
//   std::vector<double> diagonal;
//   std::vector<double> initial_vector;
// diagonal.reserve(wf.size());
// initial_vector.reserve(wf.size());
// for (const auto& term : wf.get_terms()) {
//   const auto& det = term.det;
//   diagonal.push_back(hamiltonian(det, det));
//   initial_vector.push_back(term.coef);
// }
// HelperStrings helper_strings(wf.get_dets());
// std::function<std::vector<double>(std::vector<double>)> apply_hamiltonian_func =
//     std::bind(&Solver::apply_hamiltonian, this, std::placeholders::_1, helper_strings);

// Davidson davidson(diagonal, apply_hamiltonian_func, wf.size());
// if (Parallel::get_id() == 0) davidson.set_verbose(true);
// davidson.diagonalize(initial_vector);

// double energy_var = davidson.get_lowest_eigenvalue();
// const auto& coefs_new = davidson.get_lowest_eigenvector();
// wf.set_coefs(coefs_new);
// wf.sort_by_coefs();

// return energy_var;
// }