#include "spin_det.h"

void SpinDet::set_orb(const int orb_id, const bool occ) {
  if (elecs.empty()) {
    elecs.push_back(orb_id);
  } else if (orb_id < elecs.front()) {
    elecs.insert(elecs.begin(), orb_id);
  } else if (orb_id > elecs.back()) {
    elecs.insert(elecs.end(), orb_id);
  } else {
    auto it = std::lower_bound(elecs.begin(), elecs.end(), orb_id);
    const Orbital orb_id_cast = static_cast<Orbital>(orb_id);
    if (occ) {
      if (*it == orb_id_cast) return;  // Already occupied.
      elecs.insert(it, orb_id_cast);
    } else {
      if (*it != orb_id_cast) return;
      elecs.erase(it);
    }
  }
}

void SpinDet::from_eor(const SpinDet& lhs, const SpinDet& rhs) {
  // Find the orbitals where lhs and rhs differ from each other.
  // Store in ascending order.
  const auto& lhs_elecs = lhs.elecs;
  const auto& rhs_elecs = rhs.elecs;
  const int lhs_size = lhs_elecs.size();
  const int rhs_size = rhs_elecs.size();
  int lhs_ptr = 0;
  int rhs_ptr = 0;
  elecs.clear();
  elecs.reserve(lhs_size + rhs_size);
  while (lhs_ptr < lhs_size && rhs_ptr < rhs_size) {
    if (lhs_elecs[lhs_ptr] < rhs_elecs[rhs_ptr]) {
      elecs.push_back(lhs_elecs[lhs_ptr]);
      lhs_ptr++;
    } else if (lhs_elecs[lhs_ptr] > rhs_elecs[rhs_ptr]) {
      elecs.push_back(rhs_elecs[rhs_ptr]);
      rhs_ptr++;
    } else {
      lhs_ptr++;
      rhs_ptr++;
    }
  }
  while (lhs_ptr < lhs_size) {
    elecs.push_back(lhs_elecs[lhs_ptr]);
    lhs_ptr++;
  }
  while (rhs_ptr < rhs_size) {
    elecs.push_back(rhs_elecs[rhs_ptr]);
    rhs_ptr++;
  }
}

bool operator==(const SpinDet& lhs, const SpinDet& rhs) { return lhs.elecs == rhs.elecs; }

std::ostream& operator<<(std::ostream& os, const SpinDet& spin_det) {
  for (const Orbital orbital : spin_det.elecs) {
    os << orbital << " ";
  }
  return os;
}
