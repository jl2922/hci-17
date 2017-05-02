#ifndef HCI_HEG_SOLVER_H_
#define HCI_HEG_SOLVER_H_

#include <boost/functional/hash.hpp>
#include "../std.h"

#include "../det/det.h"
#include "../solver/solver.h"
#include "../types.h"

class HEGSolver : public Solver {
 private:
  double rcut_var;
  std::size_t n_orbs_var;
  double rcut_pt;
  std::size_t n_orbs_pt;
  std::vector<std::size_t> n_orbs_pts;  // Corresponding to rcut_pts.
  double k_unit;
  double H_unit;
  std::vector<Int3> k_points;  // O(k_points).
  std::unordered_map<Int3, std::size_t, boost::hash<Int3>> k_lut;  // O(k_points).
  std::unordered_map<TinyInt3, std::vector<TinyInt3Double>, boost::hash<TinyInt3>>
      same_spin_hci_queue;  // O(k_points^2).
  std::vector<TinyInt3Double> opposite_spin_hci_queue;  // O(k_points).

  static HEGSolver get_instance() {
    static HEGSolver heg_solver;
    return heg_solver;
  }

  void solve() override;

  void setup() override;

  void generate_k_points(const double rcut);

  void generate_hci_queue(const double rcut);

  void save_variation_result();

  bool load_variation_result();

  PTCategory get_category(const Det&, const double);

  void perturbation();

  double hamiltonian(const Det&, const Det&) const override;

  std::list<Det> find_connected_dets(const Det&, const double eps) const override;

 public:
  static void run() { HEGSolver::get_instance().solve(); }
};

#endif