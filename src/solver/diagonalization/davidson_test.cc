#include "davidson.h"
#include "gtest/gtest.h"

// Test with a Hilbert matrix.
class HilbertSystem {
 public:
  HilbertSystem(int n) { this->n = n; }

  double get_hamiltonian(int i, int j) {
    const double GAMMA = 10.0;
    if (i == j) return -1.0 / (2 * i + 1);
    return -1.0 / GAMMA / (i + j + 1);
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
};

TEST(DavidsonTest, HilbertSystem) {
  const int N = 1000;
  HilbertSystem hamiltonian(N);
  std::vector<double> diagonal(N);
  for (std::size_t i = 0; i < N; i++) diagonal[i] = hamiltonian.get_hamiltonian(i, i);

  std::function<std::vector<double>(std::vector<double>)> apply_hamiltonian =
      std::bind(&HilbertSystem::apply_hamiltonian, &hamiltonian, std::placeholders::_1);

  Davidson davidson(diagonal, apply_hamiltonian, N);

  const std::vector<double> expected_eigenvalues(
      {-1.00956719, -0.3518051, -0.23097854, -0.17336724, -0.13218651});
  const std::vector<std::vector<double>> expected_eigenvectors(
      {{0.99292536, 0.08026708, 0.04720676, 0.03412438, 0.02694173},
       {0.10953429, -0.90014126, -0.22310872, -0.14701356, -0.11467862},
       {0.04208261, 0.42251014, -0.54880665, -0.25894711, -0.19182954},
       {0.00259482, -0.02195869, -0.78985725, 0.1487066, 0.10266289},
       {0.01203533, 0.04023094, 0.09953056, -0.90203616, -0.06584302}});

  // Check eigenvalue and eigenvector with reference values from exact diagonalization.
  std::vector<double> initial_vector(N, 0.0);
  initial_vector[0] = 1.0;
  davidson.diagonalize(initial_vector);
  const double lowest_eigenvalue = davidson.get_lowest_eigenvalue();
  EXPECT_NEAR(lowest_eigenvalue, expected_eigenvalues[0], 1.0e-6);
  const std::vector<double> lowest_eigenvector = davidson.get_lowest_eigenvector();
  for (int i = 0; i < 5; i++) {
    EXPECT_NEAR(lowest_eigenvector[i], expected_eigenvectors[0][i], 1.0e-4);
  }

  // Check
  std::vector<std::vector<double>> initial_vectors(5);
  for (int i = 0; i < 5; i++) {
    initial_vectors[i].assign(N, 0.0);
    initial_vectors[i][i] = 1.0;
  }
  davidson.diagonalize(initial_vectors);
  davidson.set_verbose(true);
  const std::vector<double> lowest_eigenvalues = davidson.get_lowest_eigenvalues();
  for (int i = 0; i < 5; i++) {
    EXPECT_NEAR(lowest_eigenvalues[i], expected_eigenvalues[i], 1.0e-6);
  }
}