#include "heg_solver.h"

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
      Time::end(eps_var_event);
      exit(1);
    }
    Time::end(rcut_var_event);
  }
  Time::end("variation all");

  Time::start("perturbation all");
  Time::end("perturbation all");
}