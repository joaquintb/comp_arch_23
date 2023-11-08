#include <simulation.hpp>

double const Simulation::get_sm_len() {
    return this->sm_len;
};

double const Simulation::get_mass() {
    return this->mass;
};

double const Simulation::get_num_p() {
    return this->num_p;
};

double const Simulation::get_ppm() {
    return this->ppm;
};


Simulation::Simulation(double ppm, double num_p) {
    this->ppm = ppm;
    this->num_p = num_p;
    
    this->sm_len = this->radius / this->ppm;
    this->mass = this->fluid_density / std::pow(this->ppm, 3);

};
