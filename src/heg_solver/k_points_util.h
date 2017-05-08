#ifndef HCI_K_POINTS_UTIL_H_
#define HCI_K_POINTS_UTIL_H_

#include <boost/functional/hash.hpp>
#include "../std.h"

#include "../array_math.h"
#include "../types.h"

class KPointsUtil {
 public:
  static std::size_t get_n_k_points(const double rcut);
  static std::vector<Int3> generate_k_points(const double rcut);
  static std::unordered_map<Int3, std::size_t, boost::hash<Int3>> generate_k_lut(
      const std::vector<Int3>& k_points);
  static std::vector<TinyInt3> get_k_diffs(const std::vector<Int3>& k_points);
};

#endif