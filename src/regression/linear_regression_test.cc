#include "linear_regression.h"

int main() {
  std::vector<std::vector<double>> X({{5.00E-06, 2.50E-01},
                                      {2.00E-06, 4.00E-02},
                                      {1.00E-06, 1.00E-02},
                                      {5.00E-07, 2.50E-03},
                                      {2.00E-07, 4.00E-04},
                                      {1.00E-07, 1.00E-04}});
  std::vector<double> y(
      {-0.593078303, -0.592875163, -0.592828326, -0.592809134, -0.592801929, -0.592801414});

  LinearRegression lr(X, y);
  const auto& estimate = lr.get_estimate();
  const auto& stdev = lr.get_stdev();
  const auto& prob_t = lr.get_prob_t();
  for (int i = 0; i < 3; i++) {
    printf("%20.10g %20.10g %20.10g\n", estimate[i], stdev[i], prob_t[i]);
  }
  return 0;
}