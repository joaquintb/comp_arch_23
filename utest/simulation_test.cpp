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
}

TEST(SimulationUnitTest, CheckSimulationInitializationValuesTest) {
    const int ppm = 1000.0;
    const int num_p = 1000;
    const Simulation sim = Simulation(ppm, num_p);
    
    EXPECT_EQ(sim.ppm, ppm);
    EXPECT_EQ(sim.num_p, num_p);

    EXPECT_EQ(sim.sm_len, sim.radius / sim.ppm);
    EXPECT_EQ(sim.mass, sim.fluid_density / std::pow(sim.ppm, 3));
    EXPECT_EQ(sim.n_x, floor((sim.b_max[0] - sim.b_min[0]) / sim.sm_len));
    EXPECT_EQ(sim.n_y, floor((sim.b_max[1] - sim.b_min[1]) / sim.sm_len));
    EXPECT_EQ(sim.n_z, floor((sim.b_max[2] - sim.b_min[2]) / sim.sm_len));
    EXPECT_EQ(sim.num_blocks, sim.n_x * sim.n_y * sim.n_z);
    EXPECT_EQ(sim.sm_len_sq, std::pow(sim.sm_len, 2));
    EXPECT_EQ(sim.fact_5_acc, Simulation::mew * sim.mass);

    EXPECT_EQ(sim.size_blocks[0], (sim.b_max[0] - sim.b_min[0]) / sim.n_x);
    EXPECT_EQ(sim.size_blocks[1], (sim.b_max[1] - sim.b_min[1]) / sim.n_y);
    EXPECT_EQ(sim.size_blocks[2], (sim.b_max[2] - sim.b_min[2]) / sim.n_z);

    EXPECT_EQ(sim.common_factor_acc, Simulation::forty_five_over_pi / std::pow(sim.sm_len, 6.0));
    EXPECT_EQ(sim.fact_1_acc, sim.mass * sim.p_s * 0.5);
    EXPECT_EQ(sim.sm_len_six, std::pow(sim.sm_len, 6.0));
    EXPECT_EQ(sim.prefactor_dens, (315.0 * sim.mass) / (64.0 * std::numbers::pi * pow(sim.sm_len, 9.0)));
}

TEST(SimulationUnitTest, CheckPositiveParticlesTest) {
    const int ppm = 1000.0;
    const int num_p = 1000;
    const Simulation sim = Simulation(ppm, num_p);

    EXPECT_EQ(sim.check_positive_particles(), true);
}
