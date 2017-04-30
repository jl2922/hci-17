#ifndef HCI_SOLVER_H_
#define HCI_SOLVER_H_

#include "../std.h"

#include "../det/det.h"
#include "../wavefunction/wavefunction.h"
#include "helper_strings.h"

class Solver {
 protected:
  std::size_t n_up;
  std::size_t n_dn;
  double max_abs_H;
  Wavefunction wf;
  double energy_hf;
  double energy_var;
  double energy_pt;
  double eps_var;
  double eps_pt;

  virtual void solve() {
    setup();
    variation();
    perturbation();
  }  // Controls the general solving procedure.

  virtual void setup() {}  // Generate hci queues and starting wavefunction.

  virtual double hamiltonian(const Det&, const Det&) const = 0;

  virtual std::list<Det> find_connected_dets(const Det&, const double eps) const = 0;

  std::vector<double> apply_hamiltonian(const std::vector<double>&, HelperStrings&);

  void variation();

  double diagonalize();

  virtual void perturbation() {}

  std::list<Det> find_next_dets();
};

#endif