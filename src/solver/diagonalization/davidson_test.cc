#include "davidson.h"

#include <cassert>
#include <cstddef>
#include <cstdio>
#include <functional>
#include <list>
#include <vector>

// Test with a Hilbert matrix.
class TestHamiltonian {
 public:
  TestHamiltonian(int n, int gamma) {
    this->n = n;
    this->gamma = gamma;
  }

  double get_hamiltonian(int i, int j) {
    if (i == j) return -1.0 / (2 * i + 1);
    return -1.0 / gamma / (i + j + 1);
  }

  std::vector<double> apply_hamiltonian(const std::vector<double>& v) {
    std::vector<double> Hv(n, 0.0);
    for (int i = 0; i < n; i++) {
      Hv[i] += get_hamiltonian(i, i) * v[i];
      for (int j = i + 1; j < n; j++) {
        double h_ij = get_hamiltonian(i, j);
        Hv[i] += h_ij * v[j];
        Hv[j] += h_ij * v[i];
      }
    }
    return Hv;
  }

 private:
  int n;
  int gamma;
};

int main(int argc, char** argv) {
  int N = 1000;
  const double GAMMA = 10.0;
  TestHamiltonian hamiltonian(N, GAMMA);
  printf("Hilbert Matrix (Gamma = %.2f, N = %d)\n", GAMMA, N);

  std::vector<double> diagonal(N);
  for (std::size_t i = 0; i < N; i++) diagonal[i] = hamiltonian.get_hamiltonian(i, i);

  std::function<std::vector<double>(std::vector<double>)> apply_hamiltonian =
      std::bind(&TestHamiltonian::apply_hamiltonian, &hamiltonian, std::placeholders::_1);

  Davidson davidson(diagonal, apply_hamiltonian, N);

  std::vector<double> initial_vector(N, 0.0);
  initial_vector[0] = 1.0;

  davidson.diagonalize(initial_vector);

  const double lowest_eigenvalue = davidson.get_lowest_eigenvalue();
  assert(fabs(lowest_eigenvalue - (-1.00956710)) < 1.0e-6);
  printf("Lowest Eigenvalue: %.10f\n", davidson.get_lowest_eigenvalue());

  return 0;
}