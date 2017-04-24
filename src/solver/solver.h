#ifndef HCI_SOLVER_H_
#define HCI_SOLVER_H_

class Solver {
 protected:
  void solve();
  virtual void setup() = 0;
};

#endif