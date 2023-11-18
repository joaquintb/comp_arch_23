#include "progargs.hpp"

int handle_num_args(int argc, std::string const & steps) {
  // (-1): invalid number of arguments
  if (argc - 1 != 3) {
    std::cerr << "Error: Invalid number of arguments: " << std::to_string(argc - 1) << "\n";
    std::exit(-1);
  }
  // (-1): first argument integer number
  int n_steps = 0;
  try {
    n_steps = std::stoi(steps);
  } catch (std::invalid_argument const & e) {
    std::cerr << "Error: time steps must be numeric."
              << "\n";
    std::exit(-1);
  }
  // (-2): invalid number of time steps
  if (n_steps < 0) {
    std::cerr << "Error: Invalid number of time steps."
              << "\n";
    std::exit(-2);
  }
  return n_steps;
}

void check_files(std::ifstream & inputFile, std::string const & inputName,
                 std::ofstream & outputFile, std::string const & outputName) {
  if (!inputFile.is_open()) {
    std::cerr << "Error: Cannot open " << inputName << " for reading.\n";
    std::exit(-3);
  }

  if (!outputFile.is_open()) {
    std::cerr << "Error: Cannot open " << outputName << " for writing.\n";
    std::exit(-4);
  }
}
