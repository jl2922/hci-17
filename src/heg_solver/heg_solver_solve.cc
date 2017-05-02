#include "heg_solver.h"

#include <boost/format.hpp>

#include "../config.h"
#include "../parallel.h"
#include "../time.h"

void HEGSolver::solve() {
  printf("Proc %d running on %s\n", Parallel::get_id(), Parallel::get_host().c_str());

  Time::start("variation stage");
  this->n_up = Config::get<std::size_t>("n_up");
  this->n_dn = Config::get<std::size_t>("n_dn");
  const auto& rcut_vars = Config::get_array<double>("rcut_vars");
  const auto& eps_vars = Config::get_array<double>("eps_vars");
  const auto& rcut_pts = Config::get_array<double>("rcut_pts");
  const auto& eps_pts = Config::get_array<double>("eps_pts");
  for (const double rcut_var : rcut_vars) {
    std::string rcut_var_event = str(boost::format("var with rcut_var: %#.4g") % rcut_var);
    Time::start(rcut_var_event);
    this->rcut_var = rcut_var;
    setup();
    for (const double eps_var : eps_vars) {
      std::string eps_var_event = str(boost::format("var with eps_var: %#.4g") % eps_var);
      Time::start(eps_var_event);
      this->eps_var = eps_var;
      if (!load_variation_result()) {
        variation();
        save_variation_result();
      }
      Time::end(eps_var_event);
    }
    Time::end(rcut_var_event);
  }
  Time::end("variation stage");

  Time::start("perturbation stage");
  for (const double rcut_pt : rcut_pts) {
    generate_k_points(rcut_pt);
    n_orbs_pts.push_back(k_points.size());
  }
  for (const double rcut_var : rcut_vars) {
    std::string rcut_var_event = str(boost::format("pt with rcut_var: %#.4g") % rcut_var);
    Time::start(rcut_var_event);
    this->rcut_var = rcut_var;
    for (const double eps_var : eps_vars) {
      std::string eps_var_event = str(boost::format("pt with eps_var: %#.4g") % eps_var);
      Time::start(eps_var_event);
      this->eps_var = eps_var;
      assert(load_variation_result());
      perturbation();
      Time::end(eps_var_event);
    }
    Time::end(rcut_var_event);
  }
  Time::end("perturbation stage");
}

void HEGSolver::save_variation_result() {
  if (Parallel::get_id() != 0) return;
  std::ofstream var_file;
  std::string filename = str(boost::format("var_%.5f_%.3f.txt") % eps_var % rcut_var);
  var_file.open(filename);
  var_file << boost::format("%.15g %.15g\n") % energy_hf % energy_var;
  var_file << boost::format("%d %d %d\n") % n_up % n_dn % wf.size();
  for (const auto& term : wf.get_terms()) {
    var_file << boost::format("%.15g\n") % term.coef;
    var_file << term.det.up << std::endl << term.det.dn << std::endl;
  }
  var_file.close();
  printf("Variation result saved to: %s\n", filename.c_str());
}

bool HEGSolver::load_variation_result() {
  std::ifstream var_file;
  std::string filename = str(boost::format("var_%.5f_%.3f.txt") % eps_var % rcut_var);
  std::size_t wf_size;
  int orb_id;
  double coef;
  var_file.open(filename);
  if (!var_file.is_open()) return false;
  var_file >> energy_hf >> energy_var;
  var_file >> n_up >> n_dn >> wf_size;
  wf.clear();
  for (std::size_t i = 0; i < wf_size; i++) {
    var_file >> coef;
    Det det;
    for (std::size_t j = 0; j < n_up; j++) {
      var_file >> orb_id;
      det.up.set_orb(orb_id, true);
    }
    for (std::size_t j = 0; j < n_dn; j++) {
      var_file >> orb_id;
      det.dn.set_orb(orb_id, true);
    }
    wf.append_term(det, coef);
  }
  var_file.close();
  if (Parallel::get_id() == 0) printf("Variation result loaded from: %s\n", filename.c_str());
  return true;
}