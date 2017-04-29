#include "heg_solver.h"

#include <boost/format.hpp>

#include "../config.h"
#include "../parallel.h"
#include "../time.h"

void HEGSolver::solve() {
  printf("Proc %d running on %s\n", Parallel::get_id(), Parallel::get_host().c_str());

  Time::start("variation all");
  this->n_up = Config::get<int>("n_up");
  this->n_dn = Config::get<int>("n_dn");
  const auto& rcuts_var = Config::get_array<double>("rcuts_var");
  for (const double rcut_var : rcuts_var) {
    std::string rcut_var_event = "rcut_var: " + std::to_string(rcut_var);
    Time::start(rcut_var_event);
    this->rcut_var = rcut_var;
    setup();
    const auto& epss_var = Config::get_array<double>("epss_var");
    for (const double eps_var : epss_var) {
      std::string eps_var_event = "eps_var: " + std::to_string(eps_var);
      Time::start(eps_var_event);
      this->eps_var = eps_var;
      variation();
      dump_variation_result();
      Time::end(eps_var_event);
    }
    Time::end(rcut_var_event);
  }
  Time::end("variation all");

  Time::start("perturbation all");
  Time::end("perturbation all");
}

void HEGSolver::dump_variation_result() {
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