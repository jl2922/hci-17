#ifndef HCI_SPIN_DET_H_
#define HCI_SPIN_DET_H_

#include <algorithm>
#include <cstddef>
#include <cstdint>

#include "../types.h"

class SpinDet {
  typedef uint16_t Orbital;

 private:
  vector<Orbital> elecs;

 public:
  void set_orb(const int orb_id, const bool occ);

  bool get_orb(const int orb_id) const {
    return std::binary_search(elecs.begin(), elecs.end(), orb_id);
  }

  int get_n_elecs() const { return elecs.size(); }

  void from_eor(const SpinDet&, const SpinDet&);

  vector<int> get_elec_orbs() const {
    vector<int> elecs_cast;
    elecs_cast.reserve(elecs.size());
    for (const auto elec : elecs) elecs_cast.push_back(static_cast<int>(elec));
    return elecs_cast;
  }

  void print() {
    for (std::size_t i = 0; i < elecs.size(); i++) printf("%d ", elecs[i]);
    printf("\n");
  }

  friend bool operator==(const SpinDet&, const SpinDet&);
};

bool operator==(const SpinDet&, const SpinDet&);

#endif