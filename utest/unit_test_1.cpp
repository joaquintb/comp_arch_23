#include <gtest/gtest.h>

#include "../sim/block.hpp"
#include "../sim/grid.hpp"
#include "../sim/simulation.hpp"
#include "../sim/progargs.hpp"
#include <cmath>


TEST(SimulationUnitTest, CheckConstantsTest) {
    const int ppm = 1000.0;
    const int num_p = 1000;
    const Simulation sim = Simulation(ppm, num_p);

    // Check Constants 
    EXPECT_EQ(sim.radius, 1.695);
    EXPECT_EQ(sim.fluid_density, 1000.0);
    EXPECT_EQ(sim.p_s, 3.0);
    EXPECT_EQ(sim.s_c, 30000);
    EXPECT_EQ(sim.d_v, 128.0);
    EXPECT_EQ(sim.mew, 0.4);
    EXPECT_EQ(sim.d_p, 0.0002);
    EXPECT_EQ(sim.delta_t, 0.001);

    EXPECT_EQ(sim.b_max[0], 0.065);
    EXPECT_EQ(sim.b_max[1], 0.1);
    EXPECT_EQ(sim.b_max[2], 0.065);
    EXPECT_EQ(sim.b_min[0], -0.065);
    EXPECT_EQ(sim.b_min[1], -0.08);
    EXPECT_EQ(sim.b_min[2], -0.065);

    // Check smoothing length and mass 
    EXPECT_EQ(sim.sm_len, sim.radius / sim.ppm);
    EXPECT_EQ(sim.mass, sim.fluid_density / std::pow(sim.ppm, 3));
}

// TEST(SimulationUnitTest, CheckPositiveParticlesTest) {
//     const int ppm = 1000.0;
//     const int num_p = 1000;
//     const Simulation sim = Simulation(ppm, num_p);

//     EXPECT_NO_THROW({sim.check_positive_particles();});
// }
