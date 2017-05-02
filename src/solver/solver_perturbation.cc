#include "solver.h"

#include "../parallel.h"

unsigned long long Solver::estimate_n_pt_dets(const double eps_pt) {
  const std::size_t SAMPLE_INTERVAL_MIN = 100;
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