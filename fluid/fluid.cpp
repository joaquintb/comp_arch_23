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

int main (int argc, char **argv) {
    // Grid initialization 
    // Read binary file 
    std::string filename = "small.fld";

    std::ifstream inputFile("../../small.fld", std::ios::binary);

    if (!inputFile.is_open()) {
        std::cerr << "Unable to open the file." << std::endl;
        return 1;
    }
    float ppm; 
    int num_p;

    inputFile.read(reinterpret_cast<char*>(&ppm), sizeof(float));
    inputFile.read(reinterpret_cast<char*>(&num_p), sizeof(int));

    // Create simulation 
    Simulation sim = Simulation(ppm, num_p);

    int const nx = floor((sim.b_max[0] - sim.b_min[0]) / sim.get_sm_len());
    int const ny = floor((sim.b_max[1] - sim.b_min[1]) / sim.get_sm_len());
    int const nz = floor((sim.b_max[2] - sim.b_min[2]) / sim.get_sm_len());

    // Create grid
    Grid grid = Grid(nx, ny, nz);

    // Add blocks to grid 
    grid.populate(sim, inputFile);

    return 0; 
}