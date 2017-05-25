#ifndef HCI_SPIN_DET_H_
#define HCI_SPIN_DET_H_

#include "../std.h"

#include "../types.h"

class SpinDet {
 public:
  enum EncodeScheme { FIXED, VARIABLE };

  void set_orb(const int orb_id, const bool occ);

  bool get_orb(const int orb_id) const {
    return std::binary_search(elecs.begin(), elecs.end(), orb_id);
  }

  std::size_t get_n_elecs() const { return elecs.size(); }

  void from_eor(const SpinDet&, const SpinDet&);

  const Orbitals& get_elec_orbs() const { return elecs; }

  const Orbitals encode(const EncodeScheme scheme = VARIABLE) const {
    if (scheme == FIXED) return elecs;
    return encodeVariable();
  }

  void decode(const Orbitals& code, const EncodeScheme scheme = VARIABLE) {
    if (scheme == FIXED)
      elecs = code;
    else
      decodeVariable(code);
  }

  friend bool operator==(const SpinDet&, const SpinDet&);
  friend std::ostream& operator<<(std::ostream&, const SpinDet&);

 private:
  Orbitals elecs;

  const Orbitals encodeVariable() const;

  void decodeVariable(const Orbitals& code);
};

bool operator==(const SpinDet&, const SpinDet&);
std::ostream& operator<<(std::ostream&, const SpinDet&);

#endif