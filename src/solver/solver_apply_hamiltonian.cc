#include "solver.h"

#include "../parallel.h"
#include "helper_strings.h"

std::vector<double> Solver::apply_hamiltonian(
    const std::vector<double>& vec, HelperStrings& helper_strings) {
  std::size_t n = vec.size();
  assert(n == wf.size());
  std::vector<double> res(n, 0.0);
  const auto& dets = wf.get_dets();
  for (std::size_t i = 0; i < n; i++) {
    if (i % Parallel::get_n() != static_cast<std::size_t>(Parallel::get_id())) continue;
    const Det& det_i = dets[i];
    const auto& connections = helper_strings.find_potential_connections(i);
    for (std::size_t j : connections) {
      if (j < i) continue;
      const Det& det_j = dets[j];
      const double H_ij = hamiltonian(det_i, det_j);
      res[i] += H_ij * vec[j];
      if (i != j) res[j] += H_ij * vec[i];
    }
  }
  for (std::size_t i = 0; i < n; i++) Parallel::all_reduce(res[i]);
  return res;
}