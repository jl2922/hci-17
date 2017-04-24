#include "solver.h"

#include "../config.h"
#include "../parallel.h"
#include "../time.h"

void Solver::solve() {
  printf("Proc %d running on %s\n", Parallel::get_id(), Parallel::get_proc_name().c_str());

  Time::start("solve");
  const auto& rcuts_var = Config::get_array<double>("rcuts_var");
  for (const double rcut_var : rcuts_var) {
    std::string rcut_var_event = "rcut_var: " + std::to_string(rcut_var);
    Time::start(rcut_var_event);
    const auto& epss_var = Config::get_array<double>("epss_var");
    for (const double eps_var : epss_var) {
      std::string eps_var_event = "  eps_var: " + std::to_string(eps_var);
      Time::start(eps_var_event);
      // Solve PTs.
      Time::end(eps_var_event);
    }
    Time::end(rcut_var_event);
  }
  Time::end("solve");
}