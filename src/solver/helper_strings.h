#ifndef HCI_HELPER_STRINGS_H_
#define HCI_HELPER_STRINGS_H_

#include <boost/functional/hash.hpp>
#include "../std.h"

#include "../det/det.h"
#include "../det/spin_det.h"
#include "../types.h"

class HelperStrings {
 public:
  HelperStrings(const std::vector<Det>& dets) : dets(dets) {
    setup_ab();
    setup_ab_m1();
    connected.assign(dets.size(), false);
    one_up.assign(dets.size(), false);
  }

  UnsignedInts find_potential_connections(const std::size_t i);

 private:
  // alpha and beta strings, O(n_dets).
  std::unordered_map<Orbitals, std::pair<UnsignedInts, UnsignedInts>, boost::hash<Orbitals>> ab;

  // alpha-m1 and beta-m1 strings, O(n_dets * n_elecs).
  std::unordered_map<Orbitals, std::pair<UnsignedInts, UnsignedInts>, boost::hash<Orbitals>> ab_m1;

  // Variational determinants.
  const std::vector<Det> dets;

  // Whether has been included in the potential connections.
  std::vector<bool> connected;

  // Whether the variational dets are one-up excitations of the det passed in.
  std::vector<bool> one_up;

  // Setup alpha and beta.
  void setup_ab();

  // Setup alpha-m1 and beta-m1.
  void setup_ab_m1();

  // Remove helper strings that only has one index.
  void shrink(
      std::unordered_map<Orbitals, std::pair<UnsignedInts, UnsignedInts>, boost::hash<Orbitals>>&);
};

#endif