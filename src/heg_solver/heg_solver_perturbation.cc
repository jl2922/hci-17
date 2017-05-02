#include "heg_solver.h"

#include "../big_unordered_map/big_unordered_map.h"
#include "../config.h"
#include "../parallel.h"
#include "../time.h"

void HEGSolver::perturbation() {
  // Perform perturbation with smallest eps and largest rcut.
  const double rcut_pt_max = Config::get_array<double>("rcut_pts").back();
  const double eps_pt_min = Config::get_array<double>("eps_pts").back();
  generate_k_points(rcut_pt_max);
  generate_hci_queue(rcut_pt_max);
  const std::size_t n_orbs_pt = k_points.size();
  if (Parallel::get_id() == 0) {
    printf("PT with rcut_pt_max %#.4g, eps_pt_min %#.4g\n", rcut_pt_max, eps_pt_min);
    printf("Number of max perturbation orbitals: %d\n", static_cast<int>(n_orbs_pt * 2));
  }

  // Cache variation determinants.
  var_dets_set.clear();
  for (const auto& det : wf.get_dets()) var_dets_set.insert(det.encode());
  var_dets_set.rehash(var_dets_set.size() * 2);  // <20% conflict rate.

  // Setup hash table.
  // First estimate total number of pt dets through sampling;
  // Then reserve twice amount of hash buckets.
  Time::start("setup hash table");
  unsigned long long n_pt_dets_estimate = estimate_n_pt_dets(eps_pt_min);
  if (Parallel::get_id() == 0) printf("Estimated PT dets: %'llu\n", n_pt_dets_estimate);
  std::pair<PTKey, double> skeleton;  // For reducing the amount of MPI data transfer.
  std::get<0>(skeleton.first) = wf.get_terms().front().det.encode();
  BigUnorderedMap<PTKey, double, boost::hash<PTKey>> pt_sums(skeleton);
  pt_sums.reserve(n_pt_dets_estimate * 2);
  unsigned long long hash_buckets = pt_sums.bucket_count();
  if (Parallel::get_id() == 0) printf("Reserved %'llu total hash buckets.\n", hash_buckets);
  Time::end("setup hash table");

  // Search for perturbation determinants.
  Time::start("search pt dets");
  int progress = 1;  // For print.
  std::size_t i = 0;
  for (const auto& term : wf.get_terms()) {
    if ((i++) % Parallel::get_n() != static_cast<std::size_t>(Parallel::get_id())) continue;
    const auto& connected_dets = find_connected_dets(term.det, eps_pt_min / fabs(term.coef));
    for (const auto& det_a : connected_dets) {
      if (var_dets_set.count(det_a.encode()) == 1) continue;
      const double H_ai = hamiltonian(term.det, det_a);
      if (fabs(H_ai) < DBL_EPSILON) continue;
      const double partial_sum = H_ai * term.coef;
      PTCategory category = get_category(det_a, fabs(partial_sum));
      PTKey ptKey(det_a.encode(), category);
      pt_sums.async_inc(ptKey, partial_sum);
    }
    if ((i + 1) * 100 >= wf.size() * progress && Parallel::get_id() == 0) {
      const auto& local_map = pt_sums.get_local_map();
      Time::checkpoint("search pt dets");
      printf(
          "MASTER: Progress: %d%%. Local PT keys: %'lu, hash load: %.2f\n",
          progress,
          local_map.size(),
          local_map.load_factor());
      progress *= 2;
    }
  }
  pt_sums.complete_async_incs();
  unsigned long long n_pt_keys = pt_sums.size();
  if (Parallel::get_id() == 0) printf("Total PT keys: %'llu\n", n_pt_keys);
  Time::end("search pt dets");
}

PTCategory HEGSolver::get_category(const Det& det, const double partial_sum) {
  PTCategory category = 0;
  const auto& eps_pts = Config::get_array<double>("eps_pts");
  for (const double eps_pt : eps_pts) {
    if (partial_sum >= eps_pt) category++;
  }
  category <<= 2;
  const std::size_t highest_orb_up = det.up.get_elec_orbs().back();
  const std::size_t highest_orb_dn = det.dn.get_elec_orbs().back();
  for (const std::size_t n_orbs : n_orbs_pts) {
    if (highest_orb_up <= n_orbs && highest_orb_dn <= n_orbs) category++;
  }
  return category;
}