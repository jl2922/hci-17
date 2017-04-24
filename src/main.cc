#include <cstdio>
#include <stdexcept>
#include <string>

#include "config.h"
#include "heg_solver/heg_solver.h"
#include "chemistry_solver/chemistry_solver.h"
#include "parallel.h"
#include "time.h"

int main() {
  if (Parallel::get_id() == 0) {
    printf("Heat-Bath Configuration Interaction Solver\n");
  }

  Time::init();

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