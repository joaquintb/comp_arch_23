#include "../sim/block.hpp"
#include "../sim/grid.hpp"
#include "../sim/progargs.hpp"

int main(int argc, char ** argv) {
  // ------------- Input handling -------------
  std::span const args_view{argv, static_cast<std::size_t>(argc)};
  std::vector<std::string> const arguments{args_view.begin() + 1, args_view.end()};
  const int n_steps = handle_num_args(argc, arguments[0]);
  std::ifstream inputFile(arguments[1], std::ios::binary);
  std::ofstream outputFile(arguments[2], std::ios::binary | std::ios::trunc);
  check_files(inputFile, arguments[1], outputFile, arguments[2]);
  // ------------- Simulation and grid set-up -------------
  auto ppm = read_binary_value<float>(inputFile);
  auto num_p = read_binary_value<int>(inputFile);
  Simulation sim = Simulation(ppm, num_p);
  sim.check_positive_particles();  // Check number of particles > 0
  Grid grid = Grid(sim.n_x, sim.n_y, sim.n_z);
  grid.populate(sim, inputFile);
  sim.print_sim_values();  // Print initial attributes of simulation (successful input)
  grid.set_neighbors();
  inputFile.close();
  grid.simulate(n_steps, sim);
  // ------------- Output -------------
  write_binary_value(static_cast<float>(ppm), outputFile);
  write_binary_value(static_cast<int>(num_p), outputFile);
  grid.gen_output(outputFile);
  outputFile.close();
  return 0;
}
