#ifndef COMP_ARCH_23_SIMULATION_HPP
#define COMP_ARCH_23_SIMULATION_HPP

#include <iostream>
#include <vector>
#include <cmath>
#include <array>

struct Simulation{
    // Input-dependent physical params
    double ppm, sm_len, mass;
    int num_p;

    // Grid info
    int n_x, n_y, n_z;
    int num_blocks;
    std::vector<double> size_blocks;

    // Computational values shortcut
    double sm_len_sq, common_factor_acc, fact_1_acc, fact_5_acc, sm_len_six, prefactor_dens;

    // Physical params
    static constexpr double radius = 1.695;
    static constexpr double fluid_density = 1000.0;
    static constexpr double p_s = 3.0;
    static constexpr double s_c = 30000;
    static constexpr double d_v = 128.0;
    static constexpr double mew = 0.4;
    static constexpr double d_p = 0.0002;
    static constexpr double delta_t = 0.001;
    static constexpr double gravity = -9.8;
    static constexpr double delta_coll_max = 1e-10;
    static constexpr double forty_five_over_pi = 45.0 / std::numbers::pi;
    static constexpr int particle_error_code = -5;
    static constexpr std::array<double,3> b_max{0.065, 0.1, 0.065};
    static constexpr std::array<double,3> b_min{-0.065, -0.08, -0.065};

    Simulation(float ppm, int num_p);
    void check_positive_particles() const;
    void print_sim_values();
};

#endif //COMP_ARCH_23_SIMULATION_HPP
