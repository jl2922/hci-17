#include "linear_regression.h"

#include <Eigen/Dense>
#include <boost/math/special_functions/beta.hpp>
#include <stdexcept>

double get_prob_gt_t(const double t, const int df) {
  if (df <= 0) throw std::invalid_argument("non-positive degrees of freedom");
  const double x = df / (df + t * t);
  return boost::math::ibeta(0.5 * df, 0.5, x);
}

void LinearRegression::solve(
    const std::vector<std::vector<double>>& X_in, const std::vector<double>& y_in) {
  Eigen::MatrixXd X = Eigen::MatrixXd::Zero(n, p + 1);  // 1 constant.
  Eigen::VectorXd y = Eigen::VectorXd::Zero(n);
  for (std::size_t i = 0; i < n; i++) {
    y[i] = y_in[i];
    for (std::size_t j = 0; j < p; j++) X(i, j) = X_in[i][j];
    X(i, p) = 1.0;
  }
  const Eigen::MatrixXd XT = X.transpose();
  const Eigen::MatrixXd XTXInv = (XT * X).inverse();
  const Eigen::VectorXd beta = XTXInv * XT * y;
  const Eigen::VectorXd yReg = X * beta;
  const Eigen::VectorXd yErr = y - yReg;
  const double varErr = yErr.squaredNorm() / (n - p - 1);
  const Eigen::VectorXd varBeta = varErr * XTXInv.diagonal();
  estimate.resize(p + 1);
  stdev.resize(p + 1);
  t.resize(p + 1);
  prob_t.resize(p + 1);
  for (std::size_t i = 0; i <= p; i++) {
    estimate[i] = beta[i];
    stdev[i] = sqrt(varBeta[i]);
    t[i] = beta[i] / stdev[i];
    prob_t[i] = get_prob_gt_t(t[i], n - p - 1);
  }
}