#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
#include <memory>
#include <algorithm>
#include <numbers>
#include <span>

#include "../sim/block.hpp"
#include "../sim/progargs.hpp"
#include "../sim/simulation.hpp"
#include "../sim/grid.hpp"
#include "../sim/util.hpp"

int main (int argc, char **argv) {
    // Input handling
    inputTest(argc, argv);

    // Open binary file and read ppm and num_p
    std::span const args_view{argv, static_cast<std::size_t>(argc)};
    std::vector<std::string> const arguments{args_view.begin() + 1, args_view.end()};
    std::ifstream inputFile(arguments[1], std::ios::binary);

    float ppm;
    int num_p;

    inputFile.read(reinterpret_cast<char*>(&ppm), sizeof(float));
    inputFile.read(reinterpret_cast<char*>(&num_p), sizeof(int));

    // Create simulation 
    Simulation sim = Simulation(ppm, num_p);

    sim.checkValues(); // Check num_p is correct

    int const nx = floor((sim.b_max[0] - sim.b_min[0]) / sim.get_sm_len());
    int const ny = floor((sim.b_max[1] - sim.b_min[1]) / sim.get_sm_len());
    int const nz = floor((sim.b_max[2] - sim.b_min[2]) / sim.get_sm_len());
    int const n_blocks = nx * ny * nz;

    sim.printValues(n_blocks); // Print initial attributes of simulation (successful input)

    // Create grid
    Grid grid = Grid(nx, ny, nz);

    // Add blocks to grid 
    grid.populate(sim, inputFile);
    // Add neighbor blocks to all blocks
    grid.set_neighbors();
    // grid.test_neighbors();

    inputFile.close();

    // ---------------------------------------------------------------------
    grid.increase_all_dens(sim); // Precision problems of the order of 1e-28

    // Set grid to densinc trace
    std::string fileName = "../../trz/small/densinc-base-1.trz";
    std::ifstream trace(fileName, std::ios::binary);
    grid.set_to_trace(trace);
    trace.close();

    grid.trans_all_dens(sim); // Works OK!

    grid.increase_all_accs(sim); // Precision problems of the order 1e-13

    // Set grid to acctransf
    std::string fileName2 = "../../trz/small/acctransf-base-1.trz";
    std::ifstream trace2(fileName2, std::ios::binary);
    grid.set_to_trace(trace2);
    trace2.close();

    grid.part_collisions(sim); // Works OK!

    grid.motion(sim); // Precision problem order of 1e-16, 1e-17

    // Set grid to motion
    std::string fileName4 = "../../trz/small/motion-base-1.trz";
    std::ifstream trace4(fileName4, std::ios::binary);
    grid.set_to_trace(trace4);
    trace4.close();

    grid.part_box_collisions(sim); // Works OK!

    grid.repos(sim);
    grid.init_acc(); // !
    // Works OK!

    return 0;

    // ---------------------------------------------------------------------
}
