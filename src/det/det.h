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

  static Orbital get_n_orbs_used(
      const OrbitalsPair& pair,
      const SpinDet::EncodeScheme scheme = SpinDet::EncodeScheme::VARIABLE) {
    if (scheme != SpinDet::EncodeScheme::VARIABLE) {
      throw std::invalid_argument("Only VARIABLE encode scheme implemented.");
    }
    Orbital n_orbs_used_up = pair.first.size() == 1 ? pair.first.size() : pair.first.back() + 1;
    Orbital n_orbs_used_dn = pair.second.size() == 1 ? pair.second.size() : pair.second.back() + 1;
    return std::max(n_orbs_used_up, n_orbs_used_dn) * 2;
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

  std::pair<Orbitals, Orbitals> encode(
      const SpinDet::EncodeScheme scheme = SpinDet::EncodeScheme::VARIABLE) const {
    return std::pair<Orbitals, Orbitals>(up.encode(scheme), dn.encode(scheme));
  }

  void decode(
      const std::pair<Orbitals, Orbitals>& code,
      const SpinDet::EncodeScheme scheme = SpinDet::EncodeScheme::VARIABLE) {
    up.decode(code.first, scheme);
    dn.decode(code.second, scheme);
  }
};

bool operator==(const Det&, const Det&);

#endif