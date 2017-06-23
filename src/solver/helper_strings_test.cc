#include "helper_strings.h"
#include "gtest/gtest.h"

#include "../parallel.h"

TEST(HelperStringsTest, SetupAndFindConnections) {
  std::vector<Det> dets;
  Det det1, det2, det3;
  det1.up.set_orb(1, true);
  det1.up.set_orb(2, true);
  det1.dn.set_orb(1, true);
  det1.dn.set_orb(2, true);
  det2.up.set_orb(3, true);
  det2.up.set_orb(4, true);
  det2.dn.set_orb(1, true);
  det2.dn.set_orb(2, true);
  det3.up.set_orb(3, true);
  det3.up.set_orb(4, true);
  det3.dn.set_orb(5, true);
  det3.dn.set_orb(6, true);
  dets.push_back(det1);
  dets.push_back(det2);
  dets.push_back(det3);

  boost::mpi::environment env;
  Parallel::init(env);
  HelperStrings hs(dets);
  const auto& connections = hs.find_potential_connections(0);
  EXPECT_EQ(connections.size(), 2);
  EXPECT_TRUE(std::find(connections.begin(), connections.end(), 0) != connections.end());
  EXPECT_TRUE(std::find(connections.begin(), connections.end(), 1) != connections.end());
  EXPECT_TRUE(std::find(connections.begin(), connections.end(), 2) == connections.end());
}