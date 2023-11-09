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

    sim.printValues(); // Print initial attributes of simulation (successful input)

    // Create grid
    Grid grid = Grid(nx, ny, nz);

    // Add blocks to grid 
    grid.populate(sim, inputFile);
    // Add neighbor blocks to all blocks
    grid.set_neighbors();
    // grid.test_neighbors();

    inputFile.close();

    grid.increase_all_dens(sim);

    // TRACE: check if initial population of grid was correct
    // Need to specify dir wrt to fluid executable
    const std::string fileName = "../../trz/small/densinc-base-1.trz";
    std::ifstream trace(fileName, std::ios::binary);
    grid.cmp_trace(trace);
    trace.close();

    return 0;
}
