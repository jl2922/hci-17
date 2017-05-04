#include "../std.h"

class LinearRegression {
 public:
  LinearRegression(const std::vector<std::vector<double>>& X, const std::vector<double>& y) {
    assert(X.size() == y.size());
    assert(X.size() > 0);
    n = X.size();
    p = X[0].size();
    solve(X, y);
  }
  const std::vector<double>& get_estimate() { return estimate; }
  const std::vector<double>& get_stdev() { return stdev; }
  const std::vector<double>& get_prob_t() { return prob_t; }

 protected:
  std::size_t n;
  std::size_t p;
  std::vector<double> estimate;
  std::vector<double> stdev;
  std::vector<double> t;
  std::vector<double> prob_t;

  void solve(const std::vector<std::vector<double>>&, const std::vector<double>&);
};