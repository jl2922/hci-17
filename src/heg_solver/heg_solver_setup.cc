#include "heg_solver.h"

#include "../array_math.h"
#include "../config.h"
#include "../time.h"

void HEGSolver::setup() {
  const double r_s = Config::get<double>("r_s");
  const double density = 3.0 / (4.0 * M_PI * pow(r_s, 3));
  const double cell_length = pow((n_up + n_dn) / density, 1.0 / 3);
  k_unit = 2 * M_PI / cell_length;
  H_unit = 1.0 / (M_PI * cell_length);

  generate_k_points(rcut_var);
  n_orbs_var = k_points.size();

  printf("Number of variational spin orbitals: %d\n", n_orbs_var * 2);

  Time::start("Generate HCI queue.");
  generate_hci_queue(rcut_var);
  Time::end("Generate HCI queue.");

  wf.clear();
}

// Generate k vectors in ascending order of magnitude, and then x, y, z.
void HEGSolver::generate_k_points(const double rcut) {
  // Get all valid k's.
  const int n_max = floor(rcut);
  k_points.clear();
  for (int i = -n_max; i <= n_max; i++) {
    for (int j = -n_max; j <= n_max; j++) {
      for (int k = -n_max; k <= n_max; k++) {
        if (i * i + j * j + k * k > pow(rcut, 2)) continue;
        k_points.push_back(Int3({i, j, k}));
      }
    }
  }

  // Sort.
  std::stable_sort(k_points.begin(), k_points.end(), [](const Int3& a, const Int3& b) -> bool {
    return sum(square(a)) < sum(square(b));
  });

  // Generate look-up table.
  k_lut.clear();
  for (int i = 0; i < static_cast<int>(k_points.size()); i++) k_lut[k_points[i]] = i;
}

std::vector<Int3> get_k_diffs(const std::vector<Int3>& k_points) {
  // Generate all possible differences between two different k points.
  std::unordered_set<Int3, boost::hash<Int3>> k_diffs_set;
  std::vector<Int3> k_diffs;
  const int n_orbs = static_cast<int>(k_points.size());
  for (int p = 0; p < n_orbs; p++) {
    for (int q = 0; q < n_orbs; q++) {
      if (p == q) continue;
      const auto& diff_pq = k_points[q] - k_points[p];
      if (k_diffs_set.count(diff_pq) == 1) continue;
      k_diffs.push_back(diff_pq);
      k_diffs_set.insert(diff_pq);
    }
  }
  k_diffs_set.clear();

  // Sort k_diffs into ascending order so that later sorting hci queue will be faster.
  std::stable_sort(k_diffs.begin(), k_diffs.end(), [](const Int3& a, const Int3& b) -> bool {
    return norm(a) < norm(b);
  });

  return k_diffs;
}

void HEGSolver::generate_hci_queue(const double rcut) {
  same_spin_hci_queue.clear();
  opposite_spin_hci_queue.clear();
  max_abs_H = 0.0;

  // Common dependencies.
  const auto& k_diffs = get_k_diffs(k_points);

  // Same spin.
  for (const auto& diff_pq : k_diffs) {
    for (const auto& diff_pr : k_diffs) {
      const auto& diff_sr = diff_pr + diff_pr - diff_pq;  // Momentum conservation.
      if (diff_sr == 0 || norm(diff_sr) > rcut * 2) continue;
      const auto& diff_ps = diff_pr - diff_sr;
      if (diff_ps == 0) continue;
      if (sum(square(diff_pr)) == sum(square(diff_ps))) continue;
      const double abs_H = fabs(1.0 / sum(square(diff_pr)) - 1.0 / sum(square(diff_ps)));
      if (abs_H < DBL_EPSILON) continue;
      const auto& item = TinyInt3Double(cast<TinyInt>(diff_pr), abs_H * H_unit);
      same_spin_hci_queue[cast<TinyInt>(diff_pq)].push_back(item);
    }
  }
  for (auto& kv : same_spin_hci_queue) {
    auto& items = kv.second;
    std::stable_sort(
        items.begin(), items.end(), [](const TinyInt3Double& a, const TinyInt3Double& b) -> bool {
          return a.second > b.second;
        });
    max_abs_H = std::max(max_abs_H, items.front().second);
  }

  // Opposite spin.
  for (const auto& diff_pr : k_diffs) {
    const double abs_H = 1.0 / sum(square(diff_pr));
    if (abs_H < DBL_EPSILON) continue;
    const auto& item = TinyInt3Double(cast<TinyInt>(diff_pr), abs_H * H_unit);
    opposite_spin_hci_queue.push_back(item);
  }
  std::stable_sort(
      opposite_spin_hci_queue.begin(),
      opposite_spin_hci_queue.end(),
      [](const TinyInt3Double& a, const TinyInt3Double& b) -> bool { return a.second > b.second; });
  max_abs_H = std::max(max_abs_H, opposite_spin_hci_queue.front().second);
}
