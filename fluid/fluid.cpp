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

  float ppm;
  int num_p;
  inputFile.read(reinterpret_cast<char *>(&ppm), sizeof(float));
  inputFile.read(reinterpret_cast<char *>(&num_p), sizeof(int));

    // Create simulation 
    Simulation sim = Simulation(ppm, num_p);

    sim.checkValues(); // Check num_p is correct

  int const nx       = floor((sim.b_max[0] - sim.b_min[0]) / sim.get_sm_len());
  int const ny       = floor((sim.b_max[1] - sim.b_min[1]) / sim.get_sm_len());
  int const nz       = floor((sim.b_max[2] - sim.b_min[2]) / sim.get_sm_len());
  int const n_blocks = nx * ny * nz;

    sim.printValues(n_blocks); // Print initial attributes of simulation (successful input)

    // Create grid
    Grid grid = Grid(nx, ny, nz);

    // Add blocks to grid 
    grid.populate(sim, inputFile);
    // Add neighbor blocks to all blocks
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

  outputFile.write(reinterpret_cast<char const *>(&ppm), sizeof(float));
  outputFile.write(reinterpret_cast<char const *>(&num_p), sizeof(int));
  grid.gen_output(outputFile);


  outputFile.close();

    return 0;

    // ---------------------------------------------------------------------
}
