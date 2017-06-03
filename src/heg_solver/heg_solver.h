#ifndef HEG_SOLVER_H_
#define HEG_SOLVER_H_

#include "../solver/solver.h"
#include "../std.h"
#include "k_points_util.h"

typedef std::pair<KPoint, double> KPointDouble;

class HEGSolver : public Solver {
 public:
  static void run() {
    HEGSolver heg_solver;
    heg_solver.solve();
  }

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

  // Controls the overall solving procedure.
  void solve() override;

  // Setup basic quantities, k points, and hci queue.
  void setup(const double rcut);

  void generate_hci_queue(const double rcut);

  double hamiltonian(const Det& det_pq, const Det& det_rs) const override;

  std::list<Det> find_connected_dets(const Det& det, const double eps) const override;
};

#endif