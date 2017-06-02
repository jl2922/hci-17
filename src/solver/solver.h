#ifndef SOLVER_H_
#define SOLVER_H_

#include "../std.h"

class Solver {
 protected:
  std::size_t n_up;
  std::size_t n_dn;
  double max_abs_H;

  // Controls the overall solving procedure for the system.
  virtual void solve() = 0;

  virtual void variation(const double eps_var);

  virtual void generate_hf_det();
};

#endif