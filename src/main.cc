#include "config.h"
#include "heg_solver/heg_solver.h"
#include "std.h"
#include "timer.h"

int main(int argc, char** argv) {
  printf("Heat-Bath Configuration Interaction Solver\n");

  Timer::now();

  Config::load("config.json");
  const std::string& type = Config::get<std::string>("type");
  if (type == "heg") {
    HEGSolver::run();
  } else {
    throw std::invalid_argument("System type not supported");
  }

  return 0;
}