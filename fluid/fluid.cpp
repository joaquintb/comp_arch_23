#include "../sim/block.hpp"
#include "../sim/grid.hpp"
#include "../sim/progargs.hpp"
#include <cmath>

int main(int argc, char ** argv) {
  // Input handling
  std::span const args_view{argv, static_cast<std::size_t>(argc)};
  std::vector<std::string> const arguments{args_view.begin() + 1, args_view.end()};
  const int n_steps = handle_num_args(argc, arguments[0]);
  std::ifstream inputFile(arguments[1], std::ios::binary);
  std::ofstream outputFile(arguments[2], std::ios::binary | std::ios::trunc);
  check_files(inputFile, arguments[1], outputFile, arguments[2]);

  auto ppm = read_binary_value<float>(inputFile);
  auto num_p = read_binary_value<int>(inputFile);
  
  // Create simulation
  Simulation sim = Simulation(ppm, num_p);

  sim.check_positive_particles();  // Check number of particles > 0
  sim.print_sim_values();  // Print initial attributes of simulation (successful input)

  // Create grid
  Grid grid = Grid(sim.n_x, sim.n_y, sim.n_z);
  grid.populate(sim, inputFile);
  grid.set_neighbors();

  inputFile.close();
  
  for (int time_step = 0; time_step < n_steps; ++time_step) {
    if (time_step > 0) {
      grid.repos(sim);
      grid.init_acc();
    }

    grid.increase_all_dens(sim);
    grid.trans_all_dens(sim);
    grid.increase_all_accs(sim);
    grid.part_collisions(sim);
    grid.motion(sim);
    grid.part_box_collisions(sim);
  }

  write_binary_value(static_cast<float>(ppm), outputFile);
  write_binary_value(static_cast<int>(num_p), outputFile);
  grid.gen_output(outputFile);

  outputFile.close();
  
  return 0;
}
