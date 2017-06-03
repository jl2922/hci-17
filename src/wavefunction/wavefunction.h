#ifndef WAVEFUNCTION_H_
#define WAVEFUNCTION_H_

#include "../std.h"
#include "det.h"
#include "term.h"

class Wavefunction {
 private:
  std::list<Term> terms;

 public:
  Wavefunction() {}

  std::size_t size() { return terms.size(); }

  void append_term(const Det& det, const double coef) { terms.push_back(Term(det, coef)); }

  const std::list<Term>& get_terms() const { return terms; }

  void set_coefs(const std::vector<double>& coefs) {
    std::size_t i = 0;
    for (auto it = terms.begin(); it != terms.end(); it++) it->coef = coefs[i++];
  }

  void sort_by_coefs() {
    terms.sort([](const Term& a, const Term& b) -> bool { return fabs(a.coef) > fabs(b.coef); });
  }

  void clear() { terms.clear(); }
};

#endif