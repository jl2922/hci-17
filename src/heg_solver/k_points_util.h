#ifndef K_POINTS_UTIL_H_
#define K_POINTS_UTIL_H_

#include <boost/functional/hash.hpp>
#include "../std.h"

typedef std::int8_t KPointComponent;
typedef std::array<KPointComponent, 3> KPoint;

class KPointsUtil {
 public:
  static std::size_t get_n_k_points(const double rcut);
  static std::vector<KPoint> generate_k_points(const double rcut);
  static std::unordered_map<KPoint, std::size_t, boost::hash<KPoint>> generate_k_lut(
      const std::vector<KPoint>& k_points);
  static std::vector<KPoint> get_k_diffs(const std::vector<KPoint>& k_points);
};

#endif