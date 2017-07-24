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
  std::vector<long double> res_precise(n, 0.0);
  const auto& dets = wf.get_dets();
  unsigned long long n_connections = 0;
  static unsigned long long n_connections_prev = 0;
  unsigned long long same_spin_count = 0, opposite_spin_count = 0;

  for (std::size_t i = 0; i < n; i++) {
    // if (i % Parallel::get_n() != static_cast<std::size_t>(Parallel::get_id())) continue;
    const Det& det_i = dets[i];
    auto connections = helper_strings.find_potential_connections(i);
    for (std::size_t j : connections) {
      if (j < i) continue;
      const Det& det_j = dets[j];
      const double H_ij = hamiltonian(det_i, det_j);
      if (H_ij == 0) continue;
      if (det_i.up == det_j.up || det_i.dn == det_j.dn) {
        same_spin_count++;
      } else {
        opposite_spin_count++;
      }
      res_precise[i] += H_ij * vec[j];
      if (i != j) {
        res_precise[j] += H_ij * vec[i];
        n_connections += 2;
      } else {
        n_connections++;
      }
    }
  }
  Time::checkpoint("Diagonalization", "hamiltonian applied");
#ifdef __INTEL_COMPILER
  for (std::size_t i = 0; i < n; i++) Parallel::reduce_to_sum(res_precise[i]);
#else
  Parallel::reduce_to_sum(res_precise);
#endif
  Parallel::reduce_to_sum(n_connections);
  Parallel::reduce_to_sum(same_spin_count);
  Parallel::reduce_to_sum(opposite_spin_count);
  if (Parallel::get_id() == 0 && n_connections != n_connections_prev) {
    printf("Same spin: %'llu, opposite spin: %'llu\n", same_spin_count, opposite_spin_count);
    printf("Number of connections: %'llu\n", n_connections);
    n_connections_prev = n_connections;
  }
  for (std::size_t i = 0; i < n; i++) res[i] = static_cast<double>(res_precise[i]);
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
    std::unordered_set<OrbitalsPair, boost::hash<OrbitalsPair>> new_dets_set;
    for (const auto& term : wf.get_terms()) {
      const auto& connected_dets = find_connected_dets(term.det, fabs(0.1 * eps_var / term.coef));
      for (const auto& new_det : connected_dets) {
        if (var_dets_set.count(new_det.encode()) == 0 &&
            new_dets_set.count(new_det.encode()) == 0) {
          new_dets_set.insert(new_det.encode());
          new_dets.push_back(new_det);
        }
      }
    }

    if (Parallel::get_id() == 0) {
      printf("Number of new dets: %'llu\n", static_cast<BigUnsignedInt>(new_dets.size()));
    }

    const auto& filtered_dets = filter_dets(new_dets, eps_var);
    for (const auto& filtered_det : filtered_dets) {
      var_dets_set.insert(filtered_det.encode());
      wf.append_term(filtered_det, 0.0);
    }
    if (Parallel::get_id() == 0) {
      printf("Number of filtered dets: %'llu\n", static_cast<BigUnsignedInt>(filtered_dets.size()));
      printf("Number of total dets: %'llu\n", static_cast<BigUnsignedInt>(var_dets_set.size()));
    }

    energy_var_new = diagonalize(filtered_dets.size() > 0 ? 5 : 10);
    if (Parallel::get_id() == 0) printf("Variation energy: %#.15g Ha\n", energy_var_new);
    Time::end("Variation Iteration: " + std::to_string(iteration));
    iteration++;
  }

  energy_var = energy_var_new;
  if (Parallel::get_id() == 0) printf("Final variation energy: %#.15g Ha\n", energy_var);
}

std::list<Det> Solver::filter_dets(const std::list<Det>& new_dets, const double eps_var) {
  Time::start("Filter");
  std::vector<Det> dets = wf.get_dets();
  std::vector<double> coefs = wf.get_coefs();
  HelperStrings helper_strings(dets);
  std::list<Det> filtered_dets;
  Time::checkpoint("Filter", "helper strings generated");
  for (const Det& det_i : new_dets) {
    double pt_sum = 0.0;
    auto connections = helper_strings.find_potential_connections(det_i);
    for (const auto& j : connections) {
      const Det& det_j = dets[j];
      const double H_ij = hamiltonian(det_i, det_j);
      pt_sum += H_ij * coefs[j];
    }
    Parallel::reduce_to_sum(pt_sum);
    if (fabs(pt_sum) > eps_var) filtered_dets.push_back(det_i);
  }
  Time::end("Filter");
  return filtered_dets;
}

double Solver::diagonalize(std::size_t max_iterations) {
  std::vector<double> diagonal;
  std::vector<double> initial_vector;
  diagonal.reserve(wf.size());
  initial_vector.reserve(wf.size());
  for (const auto& term : wf.get_terms()) {
    const auto& det = term.det;
    diagonal.push_back(hamiltonian(det, det));
    initial_vector.push_back(term.coef);
  }
  Time::start("Diagonalization");
  HelperStrings helper_strings(wf.get_dets());
  Time::checkpoint("Diagonalization", "helper strings generated");
  std::function<std::vector<double>(std::vector<double>)> apply_hamiltonian_func =
      std::bind(&Solver::apply_hamiltonian, this, std::placeholders::_1, helper_strings);

  Davidson davidson(diagonal, apply_hamiltonian_func, wf.size());
  if (Parallel::get_id() == 0) davidson.set_verbose(true);
  davidson.diagonalize(initial_vector, max_iterations);
  Time::end("Diagonalization");
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
