#ifndef HCI_HELPER_STRINGS_H_
#define HCI_HELPER_STRINGS_H_

#include <boost/functional/hash.hpp>
#include "../std.h"

#include "../det/det.h"
#include "../det/spin_det.h"
#include "../types.h"

class HelperStrings {
 private:
  std::unordered_map<Orbitals, std::pair<Ints, Ints>, boost::hash<Orbitals>> ab;
  std::unordered_map<Orbitals, std::pair<Ints, Ints>, boost::hash<Orbitals>> ab_m1;
  const std::vector<Det> dets;
  std::vector<bool> connected;
  std::vector<bool> one_up;

  // Setup alpha and beta.
  void setup_ab();

  // Setup alpha-m1 and beta-m1.
  void setup_ab_m1();

  // Remove helper strings that only has one index.
  void shrink(
      std::unordered_map<Orbitals, std::pair<Ints, Ints>, boost::hash<Orbitals>>& helper_strings);

 public:
  HelperStrings(const std::vector<Det>& dets) : dets(dets) {
    setup_ab();
    setup_ab_m1();
    connected.assign(dets.size(), false);
    one_up.assign(dets.size(), false);
  }

  Ints find_potential_connections(const int i);
};

#endif