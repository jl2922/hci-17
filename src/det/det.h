#ifndef HCI_DET_H_
#define HCI_DET_H_

#include "../std.h"

#include "../types.h"
#include "spin_det.h"

class Det {
 public:
  SpinDet up;
  SpinDet dn;

  bool get_orb(const int orb_id, const int dn_offset) const {
    if (orb_id < dn_offset) return up.get_orb(orb_id);
    return dn.get_orb(orb_id - dn_offset);
  }

  void set_orb(const int orb_id, const int dn_offset, const bool occ) {
    if (orb_id < dn_offset) {
      up.set_orb(orb_id, occ);
    } else {
      dn.set_orb(orb_id - dn_offset, occ);
    }
  }

  void from_eor(const Det& lhs, const Det& rhs) {
    up.from_eor(lhs.up, rhs.up);
    dn.from_eor(lhs.dn, rhs.dn);
  }

  std::pair<Orbitals, Orbitals> encode() const {
    return std::pair<Orbitals, Orbitals>(up.encode(), dn.encode());
  }

  void decode(const std::pair<Orbitals, Orbitals>& pair) {
    up.decode(pair.first);
    dn.decode(pair.second);
  }
};

bool operator==(const Det&, const Det&);

#endif