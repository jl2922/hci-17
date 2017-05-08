#ifndef HCI_SPIN_DET_H_
#define HCI_SPIN_DET_H_

#include "../std.h"

#include "../types.h"

class SpinDet {
 private:
  Orbitals elecs;

 public:
  void set_orb(const int orb_id, const bool occ);

  bool get_orb(const int orb_id) const {
    return std::binary_search(elecs.begin(), elecs.end(), orb_id);
  }

  std::size_t get_n_elecs() const { return elecs.size(); }

  void from_eor(const SpinDet&, const SpinDet&);

  const Orbitals get_elec_orbs() const { return elecs; }

  Orbital get_highest_orb() const { return elecs.back(); }

  const Orbitals encode() const {
    Orbitals orbitals(elecs);
    orbitals.shrink_to_fit();
    return orbitals;
  }

  void decode(const Orbitals& elecs) { this->elecs = elecs; }

  friend bool operator==(const SpinDet&, const SpinDet&);
  friend std::ostream& operator<<(std::ostream&, const SpinDet&);
};

bool operator==(const SpinDet&, const SpinDet&);
std::ostream& operator<<(std::ostream&, const SpinDet&);

#endif