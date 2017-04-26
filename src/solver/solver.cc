#include "solver.h"

#include "../det/det.h"

void Solver::variation() {
  // Setup HF or existing wf as initial wf and evaluate energy.
  // Loop:
  //   Find connected determinants and add to wf.
  //   Evaluate new energy.
  //   Exit loop if energy change within threshold.

  if (wf.size() == 0) {
    Det& det = wf.append_det(Det(), 1.0);
    for (int i = 0; i < n_up; i++) det.up.set_orb(i, true);
    for (int i = 0; i < n_dn; i++) det.dn.set_orb(i, true);
    energy_hf = energy_var = get_H_elem(det, det);
    printf("HF energy: %.10f\n", energy_hf);
  }
}