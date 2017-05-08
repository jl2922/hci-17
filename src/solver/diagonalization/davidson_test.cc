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

  std::vector<double> initial_vector(N, 0.0);
  initial_vector[0] = 1.0;

  davidson.diagonalize(initial_vector);

  // Check eigenvalue and eigenvector with reference values from exact diagonalization.
  const double lowest_eigenvalue = davidson.get_lowest_eigenvalue();
  EXPECT_NEAR(lowest_eigenvalue, -1.00956719, 1.0e-6);
  const std::vector<double> lowest_eigenvector = davidson.get_lowest_eigenvector();
  const std::vector<double> expected_eigenvector_head(
      {9.92925358e-01, 8.02670792e-02, 4.72067559e-02});
  for (int i = 0; i < 3; i++) {
    EXPECT_NEAR(
        lowest_eigenvector[i],
        expected_eigenvector_head[i],
        fabs(expected_eigenvector_head[i] * 0.01));
  }
}