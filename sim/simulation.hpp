#ifndef COMP_ARCH_23_SIMULATION_HPP
#define COMP_ARCH_23_SIMULATION_HPP

#include <iostream>
#include <vector>
#include <cmath>

class Simulation{
public:
    double const radius = 1.695;
    double const fluid_density = 1000.0;
    double const p_s = 3.0;
    double const s_c = 30000;
    double const d_v = 128.0;
    double const mew = 0.4;
    double const d_p = 0.0002;
    double const delta_t = 0.001;
    double n_x, n_y, n_z;
    double sm_len_sq;
    // Factors for acceleration increase
    double common_factor_acc;
    double fact_1_acc;
    double fact_5_acc;


    double const get_ppm();
    int const get_num_p();
    double const get_sm_len();
    double const get_mass();

    std::vector<double> const b_max{0.065, 0.1, 0.065};      // x_max, y_max, z_max
    std::vector<double> const b_min{-0.065, -0.08, -0.065};

    std::vector<double> size_blocks;

    Simulation(float ppm, int num_p);
    void checkValues();
    void printValues(int);

private:
    double ppm, sm_len, mass; 
    int num_p;
};

#endif //COMP_ARCH_23_SIMULATION_HPP
