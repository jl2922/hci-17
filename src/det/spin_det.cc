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
  const std::size_t lhs_size = lhs_elecs.size();
  const std::size_t rhs_size = rhs_elecs.size();
  std::size_t lhs_ptr = 0;
  std::size_t rhs_ptr = 0;
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

const Orbitals SpinDet::encodeVariable() const {
  Orbitals code;
  const std::size_t n = get_n_elecs();
  code.push_back(n);
  Orbital level = 0;
  for (const Orbital orb : elecs) {
    while (level < n && level < orb) {
      code.push_back(level);
      level++;
    }
    if (orb >= n) code.push_back(orb);
    level = orb + 1;
  }
  return code;
}

void SpinDet::decodeVariable(const Orbitals& code) {
  const std::size_t n = code[0];
  elecs.clear();
  elecs.reserve(n);
  Orbital level = 0;
  for (std::size_t i = 1; i < code.size(); i++) {
    Orbital orb = code[i];
    while (level < n && level < orb) {
      elecs.push_back(level);
      level++;
    }
    if (orb >= n) elecs.push_back(orb);
    level = orb + 1;
  }
  while (level < n) {
    elecs.push_back(level);
    level++;
  }
}

bool operator==(const SpinDet& lhs, const SpinDet& rhs) { return lhs.elecs == rhs.elecs; }

std::ostream& operator<<(std::ostream& os, const SpinDet& spin_det) {
  for (const Orbital orbital : spin_det.elecs) os << orbital << " ";
  return os;
}
