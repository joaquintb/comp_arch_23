#include "simulation.hpp"

double const Simulation::get_sm_len() {
    return this->sm_len;
};

double const Simulation::get_mass() {
    return this->mass;
};

int const Simulation::get_num_p() {
    return this->num_p;
};

double const Simulation::get_ppm() {
    return this->ppm;
};

Simulation::Simulation(float ppm, int num_p) {
    this->ppm = static_cast<double>(ppm); // Cast to double 
    this->num_p = num_p;
    
    this->sm_len = this->radius / this->ppm;
    this->mass = this->fluid_density / std::pow(this->ppm, 3);

    this->n_x = floor((this->b_max[0] - this->b_min[0]) / this->get_sm_len());
    this->n_y = floor((this->b_max[1] - this->b_min[1]) / this->get_sm_len());
    this->n_z = floor((this->b_max[2] - this->b_min[2]) / this->get_sm_len());

    this->size_blocks = {(this->b_max[0] - this->b_min[0]) / this->n_x, 
                         (this->b_max[1] - this->b_min[1]) / this->n_y,
                         (this->b_max[2] - this->b_min[2]) / this->n_z};
};
