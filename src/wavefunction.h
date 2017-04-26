#ifndef HCI_WAVEFUNCTION_H_
#define HCI_WAVEFUNCTION_H_

#include <cstddef>

#include "det/det.h"
#include "types.h"

class Wavefunction {
 private:
  std::size_t n;
  list<Det> dets;
  list<double> coefs;

 public:
  Wavefunction() { n = 0; }

  std::size_t size() { return n; }

  Det& append_det(const Det& det, const double coef) {
    dets.push_back(det);
    coefs.push_back(coef);
    n++;
    return dets.back();
  }
};

#endif