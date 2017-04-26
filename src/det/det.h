#ifndef HCI_DET_H_
#define HCI_DET_H_

#include "spin_det.h"

class Det {
 public:
  SpinDet up;
  SpinDet dn;

  void from_eor(const Det& lhs, const Det& rhs) {
    up.from_eor(lhs.up, rhs.up);
    dn.from_eor(lhs.dn, rhs.dn);
  }
};

bool operator==(const Det&, const Det&);

#endif