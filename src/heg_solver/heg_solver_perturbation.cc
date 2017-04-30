#include "heg_solver.h"

// #include <boost/format.hpp>

#include "../config.h"
// #include "../parallel.h"
// #include "../time.h"

void HEGSolver::perturbation() {
  generate_k_points(rcut_pt);
  n_orbs_pt = k_points.size();
  generate_hci_queue(rcut_pt);
}