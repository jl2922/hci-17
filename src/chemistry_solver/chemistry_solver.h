#ifndef HCI_CHEMISTRY_SOLVER_H_
#define HCI_CHEMISTRY_SOLVER_H_

#include "../det/det.h"
#include "../solver/solver.h"

class ChemistrySolver : public Solver {
 private:
  static ChemistrySolver get_instance() {
    static ChemistrySolver chemistry_solver;
    return chemistry_solver;
  }

  void solve() {}
  void setup() {}
  double get_H_elem(const Det&, const Det&) const { return 0.0; };

 public:
  static void run() { ChemistrySolver::get_instance().solve(); }
};

#endif