#ifndef SOLVER_H_
#define SOLVER_H_

#include "../std.h"
#include "../wavefunction/wavefunction.h"

typedef std::pair<std::size_t, std::size_t> PQPair;

class Solver {
 protected:
  std::size_t n_up;
  std::size_t n_dn;
  Wavefunction wf;
  double max_abs_H;  // max absolute H in the hci queues.
  double energy_hf;  // Hatree Fock energy.
  std::vector<double> energy_vars;  // Current variation energy corresponding to the wf.

  // Controls the overall solving procedure for the system.
  virtual void solve() = 0;

  // Perform a single variation run with eps_var. Start from the current wavefunction.
  // If the corresponding wavefunction file already exists, load it directly.
  virtual void variation(const double eps_var);

  // Create a determinant with the first few orbitals occupied.
  virtual Det generate_hf_det();

  // Calculate the hamiltonian element between two determinants.
  virtual double hamiltonian(const Det&, const Det&) const = 0;

  // Find all the determinants connected to the determinant passed in.
  // Connected means that the hamiltonian element between them is larger than eps.
  virtual std::list<Det> find_connected_dets(const Det& det, const double eps) const = 0;

  // List all possible pq pairs.
  std::list<PQPair> get_pq_pairs(const Det& det, const std::size_t dn_offset) const;

  // Return eigenvalues in descending order.
  std::vector<double> diagonalize(const std::size_t n_states);
};

#endif