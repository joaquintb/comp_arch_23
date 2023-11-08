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

    double const get_ppm();
    double const get_num_p();
    double const get_sm_len();
    double const get_mass();

    std::vector<double> const b_max{0.065, 0.1, 0.065};      // x_max, y_max, z_max
    std::vector<double> const b_min{-0.065, -0.08, -0.065};

    Simulation(double ppm, double num_p);

private:
    double ppm, num_p, sm_len, mass; 
};