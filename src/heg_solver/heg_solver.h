#ifndef HCI_HEG_SOLVER_H_
#define HCI_HEG_SOLVER_H_

#include "../solver/solver.h"

class HEGSolver : public Solver {
 private:
  static HEGSolver get_instance() {
    static HEGSolver heg_solver;
    return heg_solver;
  }

  void setup() {
    
  }

 public:
  static void run() { HEGSolver::get_instance().solve(); }
};

#endif