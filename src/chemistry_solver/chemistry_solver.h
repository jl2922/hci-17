#ifndef HCI_CHEMISTRY_SOLVER_H_
#define HCI_CHEMISTRY_SOLVER_H_

#include "../solver/solver.h"

class ChemistrySolver : public Solver {
 private:
  static ChemistrySolver get_instance() {
    static ChemistrySolver chemistry_solver;
    return chemistry_solver;
  }

  void setup() {}

 public:
  static void run() { ChemistrySolver::get_instance().solve(); }
};

#endif