#include "solver.h"

#include <boost/functional/hash.hpp>

#include "../det/det.h"
#include "../parallel.h"
#include "../time/time.h"
#include "diagonalization/davidson.h"
#include "helper_strings.h"

Det Solver::generate_hf_det() {
  Det det;
  for (std::size_t i = 0; i < n_up; i++) det.up.set_orb(i, true);
  for (std::size_t i = 0; i < n_dn; i++) det.dn.set_orb(i, true);
  return det;
}

std::vector<double> Solver::apply_hamiltonian(
    const std::vector<double>& vec, HelperStrings& helper_strings) {
  std::size_t n = vec.size();
  assert(n == wf.size());
  std::vector<double> res(n, 0.0);
  const auto& dets = wf.get_dets();
  for (std::size_t i = 0; i < n; i++) {
    if (i % Parallel::get_n() != static_cast<std::size_t>(Parallel::get_id())) continue;
    const Det& det_i = dets[i];
    auto connections = helper_strings.find_potential_connections(i);
    std::sort(connections.begin(), connections.end(), [](const int a, const int b) {
      return a > b;  // Small ones first to reduce round-off error.
    });
    for (std::size_t j : connections) {
      if (j < i) continue;
      const Det& det_j = dets[j];
      const double H_ij = hamiltonian(det_i, det_j);
      res[i] += H_ij * vec[j];
      if (i != j) res[j] += H_ij * vec[i];
    }
  }
  for (std::size_t i = 0; i < n; i++) Parallel::reduce_to_sum(res[i]);
  return res;
}

void Solver::variation() {
  const double THRESHOLD = 1.0e-6;

  // Setup HF or existing wf as initial wf and evaluate energy.
  if (wf.size() == 0) {
    const Det& det_hf = generate_hf_det();
    wf.append_term(det_hf, 1.0);
    energy_hf = energy_var = hamiltonian(det_hf, det_hf);
    if (Parallel::get_id() == 0) printf("HF energy: %#.15g Ha\n", energy_hf);
  }

  double energy_var_new = 0.0;  // Ensures the first iteration will run.
  var_dets_set.clear();
  for (const auto& term : wf.get_terms()) var_dets_set.insert(term.det.encode());
  int iteration = 0;  // For print.
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

unsigned long long Solver::estimate_n_pt_dets(const double eps_pt) {
  const std::size_t SAMPLE_INTERVAL_MIN = 10;
  const std::size_t n = wf.size();
  unsigned long long estimation = 0;
  auto it = wf.get_terms().begin();
  const std::size_t sample_interval = std::max(n / 500, SAMPLE_INTERVAL_MIN);
  for (std::size_t i = 0; i < n; i++) {
    const auto& term = *it++;
    if ((i % (sample_interval * Parallel::get_n())) != sample_interval * Parallel::get_id()) {
      continue;
    }
    const auto& connected_dets = find_connected_dets(term.det, eps_pt / fabs(term.coef));
    for (const auto& det : connected_dets) {
      if (var_dets_set.count(det.encode()) == 0) estimation++;
    }
  }
  estimation *= sample_interval;
  Parallel::reduce_to_sum(estimation);
  return estimation;
}