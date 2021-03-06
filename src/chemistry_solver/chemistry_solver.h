#ifndef HCI_CHEMISTRY_SOLVER_H_
#define HCI_CHEMISTRY_SOLVER_H_

#include "../std.h"

#include "../det/det.h"
#include "../solver/solver.h"

class ChemistrySolver : public Solver {
 private:
  static ChemistrySolver get_instance() {
    static ChemistrySolver chemistry_solver;
    return chemistry_solver;
  }

  double hamiltonian(const Det&, const Det&) const override { return 0.0; };

  std::list<Det> find_connected_dets(const Det&, const double) const override {
    return std::list<Det>();
  }

 public:
  static void run() { ChemistrySolver::get_instance().solve(); }
};

#endif