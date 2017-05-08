#include "linear_regression.h"
#include "gtest/gtest.h"

TEST(LinearRegressionTest, TwoIndependentVariable) {
  std::vector<std::vector<double>> X({{5.00E-06, 2.50E-01},
                                      {2.00E-06, 4.00E-02},
                                      {1.00E-06, 1.00E-02},
                                      {5.00E-07, 2.50E-03},
                                      {2.00E-07, 4.00E-04},
                                      {1.00E-07, 1.00E-04}});
  std::vector<double> y(
      {-0.5930783034, -0.592875163, -0.5928283259, -0.5928091337, -0.5928019288, -0.5928014142});

  LinearRegression lr(X, y);

  const auto& estimate = lr.get_estimate();
  const auto& stdev = lr.get_stdev();
  const auto& prob_t = lr.get_prob_t();
  EXPECT_EQ(estimate.size(), 3);
  EXPECT_EQ(stdev.size(), 3);
  EXPECT_EQ(prob_t.size(), 3);

  std::vector<double> expected_estimate({-26.8764, -0.000590375, -0.592765});
  std::vector<double> expected_stdev({2.10596, 3.9716e-5, 1.4868e-6});
  std::vector<double> expected_prob_t({0.00103798, 0.00066062, 3.4797e-17});
  for (int i = 0; i < 3; i++) {
    EXPECT_NEAR(estimate[i], expected_estimate[i], fabs(expected_estimate[i] * 0.0001));
    EXPECT_NEAR(stdev[i], expected_stdev[i], fabs(expected_stdev[i] * 0.001));
    EXPECT_NEAR(prob_t[i], expected_prob_t[i], 1.0e-6);
  }
}

TEST(LinearRegressionTest, NotEnoughObservations) {
  std::vector<std::vector<double>> X({{5.00E-06, 2.50E-01}, {1.00E-07, 1.00E-04}});
  std::vector<double> y({-0.5930783034, -0.5928014142});

  EXPECT_THROW(LinearRegression lr(X, y), std::invalid_argument);
}