#include "heg_solver.h"

#include <boost/format.hpp>

#include "../big_unordered_map.h"
#include "../parallel.h"
#include "../time/time.h"

// Put keys with the same det onto the same process.
#ifndef SERIAL
template <>
std::size_t BigUnorderedMap<PTKey, double, boost::hash<PTKey>>::get_target(const PTKey& key) {
  boost::hash<OrbitalsPair> det_code_hasher;
  return proc_map[det_code_hasher(key.first) % total_proc_buckets];
}
#endif

void HEGSolver::perturbation() {
  // Perform perturbation with smallest eps and largest rcut.
  const double rcut_pt_max = rcut_pts.back();
  const double eps_pt_min = eps_pts.back();
  generate_k_points(rcut_pt_max);
  generate_hci_queue(rcut_pt_max);
  if (Parallel::get_id() == 0) {
    printf("PT with rcut_pt_max %#.4g, eps_pt_min %#.4g\n", rcut_pt_max, eps_pt_min);
    printf("Number of max perturbation orbitals: %d\n", static_cast<int>(k_points.size() * 2));
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
  skeleton.first.first = wf.get_terms().front().det.encode();
  BigUnorderedMap<PTKey, double, boost::hash<PTKey>> pt_sums(skeleton);
  pt_sums.reserve(n_pt_dets_estimate * 2);
  unsigned long long hash_buckets = pt_sums.bucket_count();
  if (Parallel::get_id() == 0) printf("Reserved %'llu total hash buckets.\n", hash_buckets);
  Time::end("setup hash table");

  // Search for perturbation determinants.
  Time::start("search pt dets");
  int progress = 1;  // For print.
  std::size_t i = 0;
  int n_conn = 0;
  for (const auto& term : wf.get_terms()) {
    if ((i++) % Parallel::get_n() != static_cast<std::size_t>(Parallel::get_id())) continue;
    const auto& connected_dets = find_connected_dets(term.det, eps_pt_min / fabs(term.coef));
    n_conn += connected_dets.size();
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
          "MASTER progress: %d%%. Local PT keys: %'lu, hash load: %.2f\n",
          progress,
          local_map.size(),
          local_map.load_factor());
      progress *= 2;
    }
  }
  printf("%d: n_conn: %d\n", Parallel::get_id(), n_conn);
  pt_sums.complete_async_incs();
  unsigned long long n_pt_keys = pt_sums.size();
  if (Parallel::get_id() == 0) printf("Total PT keys: %'llu\n", n_pt_keys);
  Time::end("search pt dets");

  Time::start("accumulate contributions");
  const auto& local_map = pt_sums.get_local_map();
  BigUnsignedInt n_pt_dets;
  std::unordered_set<OrbitalsPair, boost::hash<OrbitalsPair>> pt_dets;
  for (const auto& kv : local_map) {
    const PTKey& key = kv.first;
    pt_dets.insert(key.first);
  }
  n_pt_dets = pt_dets.size();
  Parallel::reduce_to_sum(n_pt_dets);
  if (Parallel::get_id() == 0) printf("Number of PT dets: %'llu\n", n_pt_dets);
  Det det_a;
  for (const std::size_t n_orbs_pt : n_orbs_pts) {
    std::string n_orbs_pt_event = "accumulate for n_orbs_pt: " + std::to_string(n_orbs_pt * 2);
    Time::start(n_orbs_pt_event);
    for (const double eps_pt : eps_pts) {
      std::string eps_pt_event = str(boost::format("accumulate for eps_pt: %.4g") % eps_pt);
      Time::start(eps_pt_event);
      // Accumulate local map contributions.
      const auto& related_categories = get_related_categories(n_orbs_pt, eps_pt);
      if (Parallel::get_id() == 0)
        printf("DEBUG: related categories: %lu\n", related_categories.size());
      energy_pt = 0.0;
      BigUnsignedInt n_pt_dets = 0;
      // std::unordered_set<OrbitalsPair, boost::hash<OrbitalsPair>> pt_dets;
      pt_dets.clear();
      for (const auto& kv : local_map) {
        const auto& key = kv.first;
        pt_dets.insert(key.first);
        const PTCategory category = key.second;
        if (std::find(related_categories.begin(), related_categories.end(), category) ==
            related_categories.end()) {
          continue;
        }
        bool is_smallest = true;
        double partial_sum = 0.0;
        for (const auto related_category : related_categories) {
          const PTKey related_key(key.first, related_category);
          if (local_map.count(related_key) == 1) {
            // Only the smallest one submits the contribution.
            if (related_category < category) {
              is_smallest = false;
              break;
            }
            partial_sum += local_map.at(related_key);
          }
        }
        if (is_smallest) {
          det_a.decode(key.first);
          const double H_aa = hamiltonian(det_a, det_a);
          energy_pt += pow(partial_sum, 2) / (energy_var - H_aa);
          n_pt_dets++;
        }
      }  // local_map loop.
      Parallel::reduce_to_sum(n_pt_dets);  // DEBUG.
      Parallel::reduce_to_sum(energy_pt);
      if (Parallel::get_id() == 0) {
        printf("Number of related PT dets: %'llu\n", n_pt_dets);
        printf("n_orbs_var: %d\n", get_n_orbs(rcut_var) * 2);
        printf("eps_var: %#.4g\n", eps_var);
        printf("n_orbs_pt: %d\n", static_cast<int>(n_orbs_pt * 2));
        printf("eps_pt: %#.4g\n", eps_pt);
        printf("Perturbation energy: %#.12g Ha\n", energy_pt);
        printf("Correlation Energy: %.12g Ha\n", energy_var + energy_pt - energy_hf);
        std::vector<double> parameter_set(
            {1.0 / (get_n_orbs(rcut_var) * 2), eps_var, 1.0 / (n_orbs_pt * 2), eps_pt});
        parameter_sets.push_back(parameter_set);
        printf("Number of parameter sets: %d\n", static_cast<int>(parameter_sets.size()));
        results.push_back(energy_var + energy_pt - energy_hf);
      }
      Time::end(eps_pt_event);
    }  // eps_pts loop.
    Time::end(n_orbs_pt_event);
  }  // n_orbs_pts loop.
  Time::end("accumulate contributions");
}

std::vector<PTCategory> HEGSolver::get_related_categories(
    const std::size_t n_orbs, const double eps) {
  std::vector<PTCategory> related_categories;
  for (const double eps_pt : eps_pts) {
    if (eps > eps_pt) continue;
    for (const std::size_t n_orbs_pt : n_orbs_pts) {
      if (n_orbs < n_orbs_pt) continue;
      if (Parallel::get_id() == 0)
        printf(
            "Related category: %d (%lu, %#.4g)\n",
            get_category(n_orbs_pt, eps_pt),
            n_orbs_pt,
            eps_pt);
      related_categories.push_back(get_category(n_orbs_pt, eps_pt));
    }
  }
  std::sort(related_categories.begin(), related_categories.end());
  return related_categories;
}

PTCategory HEGSolver::get_category(const Det& det, const double eps) {
  const std::size_t highest_orb_up = det.up.get_elec_orbs().back();
  const std::size_t highest_orb_dn = det.dn.get_elec_orbs().back();
  const std::size_t highest_orb = std::max(highest_orb_up, highest_orb_dn);
  return get_category(highest_orb + 1, eps);
}

PTCategory HEGSolver::get_category(const std::size_t n_orbs, const double eps) {
  // The higher 4 bits represent eps_pt categories and the lower 4 bits represent n_orbs_pt.
  PTCategory category1 = 0;
  PTCategory category2 = 0;
  for (const double eps_pt : eps_pts) {
    if (eps >= eps_pt) category1++;
  }
  category1--;
  for (const std::size_t n_orbs_pt : n_orbs_pts) {
    if (n_orbs <= n_orbs_pt) category2++;
  }
  category2--;
  return (category1 << 4) + category2;
}

int HEGSolver::get_n_orbs(const double rcut) {
  int n_orbs = 0;
  const int n_max = floor(rcut);
  for (int i = -n_max; i <= n_max; i++) {
    for (int j = -n_max; j <= n_max; j++) {
      for (int k = -n_max; k <= n_max; k++) {
        if (i * i + j * j + k * k > pow(rcut, 2)) continue;
        n_orbs++;
      }
    }
  }
  return n_orbs;
}