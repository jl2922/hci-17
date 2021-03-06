#ifndef SERIAL
#include <boost/mpi.hpp>
#endif
#include "std.h"

#include "chemistry_solver/chemistry_solver.h"
#include "config.h"
#include "heg_solver/heg_solver.h"
#include "parallel.h"
#include "time/time.h"

int main(int argc, char** argv) {
#ifndef SERIAL
  boost::mpi::environment env(argc, argv);  // For MPI 1.1.
  Parallel::init(env);
#endif
  Time::init();
  std::setlocale(LC_NUMERIC, "");

  if (Parallel::get_id() == 0) printf("Heat-Bath Configuration Interaction Solver\n");

  Config::load("config.json");
  Time::start("HCI");

  // Solve.
  const std::string& type = Config::get<std::string>("type");
  if (type == "heg") {
    HEGSolver::run();
  } else if (type == "chemistry") {
    ChemistrySolver::run();
  } else {
    throw std::invalid_argument("System type not supported");
  }

  Time::end("HCI");

  return 0;
}
