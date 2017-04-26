#ifndef HCI_SOLVER_H_
#define HCI_SOLVER_H_

#include "../det/det.h"
#include "../wavefunction.h"

class Solver {
 protected:
  int n_up;
  int n_dn;
  double max_abs_H;
  Wavefunction wf;
  double energy_hf;
  double energy_var;
  double energy_pt;

  virtual void solve() = 0;  // Controls the general solving procedure.
  virtual void setup() = 0;  // Generate hci queues and starting wavefunction.
  virtual double get_H_elem(const Det&, const Det&) const = 0;

  void variation();
};

#endif