#include "config.h"
#include "heg_solver/heg_solver.h"
#include "std.h"
#include "timer.h"

int main(int argc, char** argv) {
  printf("Heat-Bath Configuration Interaction Solver\n");

  // Specify numerical format for print.
  std::setlocale(LC_NUMERIC, "");

  // Print the current time.
  Timer::now();

  // Read and print the configuration.
  Config::load("config.json");

  // Run the corresponding solver.
  const std::string& type = Config::get<std::string>("type");
  if (type == "heg") {
    HEGSolver::run();
  } else {
    throw std::invalid_argument("System type not supported");
  }

  return 0;
}