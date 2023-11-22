#include <gtest/gtest.h>

#include "../sim/block.hpp"
#include "../sim/grid.hpp"
#include "../sim/simulation.hpp"
#include "../sim/progargs.hpp"
#include <cmath>
#include <fstream>
#include <sys/stat.h>


TEST(GridUnitTest, CompareTotalBlocksTest) {

  std::ifstream inputFile("small.fld", std::ios::binary);
  const std::ofstream outputFile("out.fld", std::ios::binary | std::ios::trunc);
  //std::ifstream trace("../trz/small/repos-base-1.trz", std::ios::binary);
  //int totalBlocks = 1;
  const int totalBlocks = 4725;

  //trace.read(reinterpret_cast<char *>(&totalBlocks), sizeof(totalBlocks));
  // ------------- Simulation and grid set-up -------------
  auto ppm = read_binary_value<float>(inputFile);
  auto num_p = read_binary_value<int>(inputFile);
  Simulation sim = Simulation(ppm, num_p);
  Grid grid = Grid(sim.n_x, sim.n_y, sim.n_z);
  grid.populate(sim, inputFile);
  //grid.set_neighbors();
  //sim.print_sim_values();  // Print initial attributes of simulation (successful input)
  inputFile.close();
  EXPECT_EQ(grid.size, totalBlocks);
  //grid.simulate(n_steps, sim);

}

// TEST(GridUnitTest, CompareParticleWithinEachBlock){
//     std::ifstream inputFile("small.fld", std::ios::binary);
//     const std::ofstream outputFile("out.fld", std::ios::binary | std::ios::trunc);
//     std::ifstream trace("../../trz/small/densinc-base-2.trz", std::ios::binary);
//     int32_t totalBlocks = 0;
//     totalBlocks = read_binary_value<int32_t>(trace);
//     std::cout << "Total Number of BLocks is: " << totalBlocks << "\n";
//     // ------------- Simulation and grid set-up -------------
//     auto ppm = read_binary_value<float>(inputFile);
//     auto num_p = read_binary_value<int>(inputFile);
//     Simulation sim = Simulation(ppm, num_p);
//     Grid grid = Grid(sim.n_x, sim.n_y, sim.n_z);
//     grid.populate(sim, inputFile);
//     //grid.set_neighbors();
//     //sim.print_sim_values();  // Print initial attributes of simulation (successful input)
//     inputFile.close();


//     for (int i = 0; i < totalBlocks; ++i) {
//     // Read the number of particles in the block from the binary file
//         int64_t numParticles = 0;
//         //double dummyValue = 0;
//         const int attributes = 9;
//         numParticles = read_binary_value<int64_t>(trace);
//         //trace.read(reinterpret_cast<char *>(&numParticles), sizeof(numParticles));
//         for(int k = 0; k < numParticles; k++){
//             for (int j = 0; j < attributes; j++){
//                 read_binary_value<double>(trace);
//             }
//         }
//         // Check it matches the number of particles of the block in the trace
//         EXPECT_EQ(grid.blocks[i].particles.size() , numParticles);
//     }
// }
