#include <gtest/gtest.h>

#include "../sim/block.hpp"
#include "../sim/grid.hpp"
#include "../sim/simulation.hpp"
#include "../sim/progargs.hpp"
#include <cmath>
#include <fstream>
#include <iostream>
#include <string>

class SimulationFunctionalTest : public testing::Test {
    public: 
        Simulation sim;
        Grid grid; 
        SimulationFunctionalTest() : sim(1000.0, 1000), grid(1, 1, 1) {}
        const int n_steps = 2000;

    protected:
        void SetUp() override {}
};

TEST_F(SimulationFunctionalTest, SmallInputFunctionalTest) {
    std::ifstream inputFile("small.fld", std::ios::binary);
    std::ofstream outputFile("out_ftest.fld", std::ios::binary | std::ios::trunc);

    auto ppm = read_binary_value<float>(inputFile);
    auto num_p = read_binary_value<int>(inputFile);

    sim = Simulation(ppm, num_p);
    grid = Grid(sim.n_x, sim.n_y, sim.n_z);

    grid.populate(sim, inputFile);
    grid.set_neighbors();
    inputFile.close();
    grid.simulate(n_steps, sim);

    write_binary_value(static_cast<float>(ppm), outputFile);
    write_binary_value(static_cast<int>(num_p), outputFile);
    grid.gen_output(outputFile);
    outputFile.close();

    // INSERT ASSERT STATEMENT HERE TO COMPARE OUTPUT FILE TO EXPECTED OUTPUT FILE 
}


TEST_F(SimulationFunctionalTest, LargeInputFunctionalTest) {

}