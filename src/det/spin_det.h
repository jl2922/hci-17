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

  int get_n_elecs() const { return elecs.size(); }

  void from_eor(const SpinDet&, const SpinDet&);

  const Orbitals& get_elec_orbs() const { return elecs; }

  const Orbitals& encode() const { return elecs; }

  void decode(const Orbitals& elecs) { this->elecs = elecs; }

  friend bool operator==(const SpinDet&, const SpinDet&);
};

bool operator==(const SpinDet&, const SpinDet&);

#endif