#ifndef HEG_SOLVER_H_
#define HEG_SOLVER_H_

#include "../solver/solver.h"
#include "../std.h"
#include "k_points_util.h"

typedef std::pair<KPoint, double> KPointDouble;

class HEGSolver : public Solver {
 public:
  static void run() { HEGSolver::get_instance().solve(); }

 private:
  double k_unit;
  double H_unit;
  std::vector<double> rcut_vars;
  std::vector<double> eps_vars;
  std::vector<double> rcut_pts;
  std::vector<double> eps_pts;

  std::vector<KPoint> k_points;  // O(k_points).
  std::unordered_map<KPoint, std::size_t, boost::hash<KPoint>> k_lut;  // O(k_points).
  std::unordered_map<KPoint, std::vector<KPointDouble>, boost::hash<KPoint>>
      same_spin_hci_queue;  // O(k_points^2).
  std::vector<KPointDouble> opposite_spin_hci_queue;  // O(k_points).

  void solve() override;

  void setup(const double rcut);

  void generate_hci_queue(const double rcut);

  static HEGSolver get_instance() {
    static HEGSolver heg_solver;
    return heg_solver;
  }
};

#endif