#ifndef HCI_HEG_SOLVER_H_
#define HCI_HEG_SOLVER_H_

#include <boost/functional/hash.hpp>

#include "../solver/solver.h"
#include "../types.h"

class HEGSolver : public Solver {
 private:
  double rcut_var;
  int n_orbs_var;
  double k_unit;
  vector<Int3> k_points;  // O(k_points).
  unordered_map<Int3, int, boost::hash<Int3>> k_lut;  // O(k_points).
  unordered_map<TinyInt3, vector<TinyInt3Double>, boost::hash<TinyInt3>>
      same_spin_hci_queue;  // O(k_points^2).
  vector<TinyInt3Double> opposite_spin_hci_queue; // O(k_points).

  static HEGSolver get_instance() {
    static HEGSolver heg_solver;
    return heg_solver;
  }

  void solve() override;
  void setup() override;
  void generate_k_points(const double rcut);
  void generate_hci_queue(const double rcut);

 public:
  static void run() { HEGSolver::get_instance().solve(); }
};

#endif