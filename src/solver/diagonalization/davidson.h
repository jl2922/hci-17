#ifndef DAVIDSON_H_
#define DAVIDSON_H_

#include <Eigen/Dense>
#include "../../std.h"

// Translated from Adam's fortran code.
class Davidson {
 public:
  Davidson(
      std::vector<double>& diagonal,
      std::function<std::vector<double>(std::vector<double>)>& apply_hamiltonian,
      const int n)
      : diagonal(diagonal), apply_hamiltonian(apply_hamiltonian) {
    this->n = n;
    diagonalized = false;
    verbose = false;
  }

  void set_verbose(const bool verbose) { this->verbose = verbose; }

  void diagonalize(const std::vector<double>& initial_vector);

  double get_lowest_eigenvalue() {
    if (!diagonalized) throw std::runtime_error("Accessing eigenvalue before diagonalization.");
    return lowest_eigenvalue;
  }

  const std::vector<double>& get_lowest_eigenvector() {
    if (!diagonalized) throw std::runtime_error("Accessing eigenvector before diagonalization.");
    return lowest_eigenvector;
  }

 private:
  // Use functional programming to allow either direct or indirect evaluation.
  std::vector<double>& diagonal;
  std::function<std::vector<double>(std::vector<double>)>& apply_hamiltonian;

  // Length in each direction.
  int n;

  // Solutions.
  double lowest_eigenvalue;
  std::vector<double> lowest_eigenvector;
  bool diagonalized;
  bool verbose;
};

#endif