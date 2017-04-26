#ifndef HCI_WAVEFUNCTION_H_
#define HCI_WAVEFUNCTION_H_

#include "std.h"

#include "det/det.h"

class Wavefunction {
 private:
  std::list<std::pair<Det, double>> terms;

 public:
  Wavefunction() {}

  std::size_t size() { return terms.size(); }

  void append_term(const Det& det, const double coef) {
    terms.push_back(std::pair<Det, double>(det, coef));
  }

  const std::list<std::pair<Det, double>>& get_terms() const { return terms; }

  Det& get_last_det() { return terms.back().first; }

  void clear() { terms.clear(); }
};

#endif