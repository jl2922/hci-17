#include "k_points_util.h"

std::vector<Int3> KPointsUtil::generate_k_points(const double rcut) {
  std::vector<Int3> k_points;
  const int n_max = floor(rcut);
  for (int i = -n_max; i <= n_max; i++) {
    for (int j = -n_max; j <= n_max; j++) {
      for (int k = -n_max; k <= n_max; k++) {
        if (i * i + j * j + k * k > pow(rcut, 2)) continue;
        k_points.push_back(Int3({i, j, k}));
      }
    }
  }
  std::stable_sort(k_points.begin(), k_points.end(), [](const Int3& a, const Int3& b) -> bool {
    return sum(square(a)) < sum(square(b));
  });
  return k_points;
}

std::size_t KPointsUtil::get_n_k_points(const double rcut) {
  std::size_t n_k_points = 0;
  const int n_max = floor(rcut);
  for (int i = -n_max; i <= n_max; i++) {
    for (int j = -n_max; j <= n_max; j++) {
      for (int k = -n_max; k <= n_max; k++) {
        if (i * i + j * j + k * k > pow(rcut, 2)) continue;
        n_k_points++;
      }
    }
  }
  return n_k_points;
}

std::unordered_map<Int3, std::size_t, boost::hash<Int3>> KPointsUtil::generate_k_lut(
    const std::vector<Int3>& k_points) {
  std::unordered_map<Int3, std::size_t, boost::hash<Int3>> k_lut;
  for (std::size_t i = 0; i < k_points.size(); i++) k_lut[k_points[i]] = i;
  return k_lut;
}

std::vector<TinyInt3> KPointsUtil::get_k_diffs(const std::vector<Int3>& k_points) {
  // Generate all possible differences between two different k points.
  std::unordered_set<TinyInt3, boost::hash<TinyInt3>> k_diffs_set;
  std::vector<TinyInt3> k_diffs;
  const std::size_t n_orbs = k_points.size();
  for (std::size_t p = 0; p < n_orbs; p++) {
    for (std::size_t q = 0; q < n_orbs; q++) {
      if (p == q) continue;
      const auto& diff_pq = cast<TinyInt>(k_points[q] - k_points[p]);
      if (k_diffs_set.count(diff_pq) == 1) continue;
      k_diffs.push_back(diff_pq);
      k_diffs_set.insert(diff_pq);
    }
  }
  k_diffs_set.clear();

  // Sort k_diffs into ascending order so that later sorting hci queue will be faster.
  std::stable_sort(
      k_diffs.begin(), k_diffs.end(), [](const TinyInt3& a, const TinyInt3& b) -> bool {
        return norm(a) < norm(b);
      });

  return k_diffs;
}
