#ifndef HCI_SOLVER_H_
#define HCI_SOLVER_H_

class Solver {
 protected:
  int n_up;
  int n_dn;
  double max_abs_H;

  virtual void solve() = 0; // Controls the general solving procedure.
  virtual void setup() = 0; // Generate hci queues and starting wavefunction.

  void variation() {
    // Setup HF or existing wf as initial wf and evaluate energy.
    // Loop:
    //   Find connected determinants and add to wf.
    //   Evaluate new energy
    //   Exit loop if energy change within threshold.
  }
};

#endif