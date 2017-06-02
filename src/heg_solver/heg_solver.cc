#include "heg_solver.h"

#include "../array_math.h"
#include "../config.h"
#include "../std.h"
#include "../timer.h"
#include "k_points_util.h"

void HEGSolver::solve() {
  // Get configurations.
  n_up = Config::get<std::size_t>("n_up");
  n_dn = Config::get<std::size_t>("n_dn");
  rcut_vars = Config::get_array<double>("rcut_vars");
  eps_vars = Config::get_array<double>("eps_vars");
  rcut_pts = Config::get_array<double>("rcut_pts");
  eps_pts = Config::get_array<double>("eps_pts");
  std::sort(rcut_vars.begin(), rcut_vars.end());
  std::sort(rcut_pts.begin(), rcut_pts.end());
  std::sort(eps_vars.begin(), eps_vars.end());
  std::sort(eps_pts.begin(), eps_pts.end());
  std::reverse(eps_vars.begin(), eps_vars.end());
  std::reverse(eps_pts.begin(), eps_pts.end());

  // Variation.
  Timer::start("variation");
  for (const double rcut_var : rcut_vars) {
    Timer::start(str(boost::format("rcut_var: %#.4g") % rcut_var));
    Timer::start("setup");
    setup(rcut_var);
    Timer::end();
    for (const double eps_var : eps_vars) {
      Timer::start(str(boost::format("eps_var: %#.4g") % eps_var));
      variation(eps_var);
      Timer::end();
    }
    Timer::end();
  }
  Timer::end();

  // Time::start("perturbation stage");
  // n_orbs_pts.clear();
  // for (const double rcut_pt : rcut_pts) {
  //   n_orbs_pts.push_back(KPointsUtil::get_n_k_points(rcut_pt) * 2);
  // }
  // // Start from the largest PT so that it fails earlier upon insufficient memory.
  // for (const double rcut_var : rcut_vars | boost::adaptors::reversed) {
  //   std::string rcut_var_event = str(boost::format("perturbation with rcut_var: %#.4g") %
  //   rcut_var);
  //   Time::start(rcut_var_event);
  //   this->rcut_var = rcut_var;
  //   for (const double eps_var : eps_vars | boost::adaptors::reversed) {
  //     std::string eps_var_event = str(boost::format("perturbation with eps_var: %#.4g") %
  //     eps_var);
  //     Time::start(eps_var_event);
  //     this->eps_var = eps_var;
  //     assert(load_variation_result());
  //     perturbation();
  //     Time::end(eps_var_event);
  //   }
  //   Time::end(rcut_var_event);
  // }
  // Time::end("perturbation stage");

  // Time::start("extrapolation");
  // extrapolate();
  // Time::end("extrapolation");
}

void HEGSolver::setup(const double rcut) {
  const double r_s = Config::get<double>("r_s");
  const double density = 3.0 / (4.0 * M_PI * pow(r_s, 3));
  const double cell_length = pow((n_up + n_dn) / density, 1.0 / 3);
  k_unit = 2 * M_PI / cell_length;
  H_unit = 1.0 / (M_PI * cell_length);

  k_points = KPointsUtil::generate_k_points(rcut);
  k_lut = KPointsUtil::generate_k_lut(k_points);

  generate_hci_queue(rcut);
}

void HEGSolver::generate_hci_queue(const double rcut) {
  same_spin_hci_queue.clear();
  opposite_spin_hci_queue.clear();
  max_abs_H = 0.0;

  // Common dependencies.
  const auto& k_diffs = KPointsUtil::get_k_diffs(k_points);

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
      const auto& item = KPointDouble(diff_pr, abs_H * H_unit);
      same_spin_hci_queue[diff_pq].push_back(item);
    }
  }
  for (auto& kv : same_spin_hci_queue) {
    auto& items = kv.second;
    std::stable_sort(
        items.begin(), items.end(), [](const KPointDouble& a, const KPointDouble& b) -> bool {
          return a.second > b.second;
        });
    max_abs_H = std::max(max_abs_H, items.front().second);
  }

  // Opposite spin.
  for (const auto& diff_pr : k_diffs) {
    const double abs_H = 1.0 / sum(square(diff_pr));
    if (abs_H < DBL_EPSILON) continue;
    const auto& item = KPointDouble(diff_pr, abs_H * H_unit);
    opposite_spin_hci_queue.push_back(item);
  }
  std::stable_sort(
      opposite_spin_hci_queue.begin(),
      opposite_spin_hci_queue.end(),
      [](const KPointDouble& a, const KPointDouble& b) -> bool { return a.second > b.second; });
  max_abs_H = std::max(max_abs_H, opposite_spin_hci_queue.front().second);
}