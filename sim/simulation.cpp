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

    this->sm_len_sq = std::pow(this->sm_len, 2);
    this->common_factor_acc = 45.0 / (std::numbers::pi * std::pow(this->sm_len, 6));
    this->fact_1_acc        = this->mass * this->p_s * 0.5;
    this->fact_5_acc        = this->mew * this->mass;
};

void Simulation::checkValues() {
    if(this->get_num_p() <= 0){
        std::cerr << "Error: Invalid number of particles: " << this->get_num_p() << ".\n";
        std::exit(-5);
    }
    // Pending (else if) if num_p == number of particles read in file
}

void Simulation::printValues(const int n_blocks) {
    std::cout << "Number of particles: " << this->get_num_p() << '\n';
    std::cout << "Particles per meter: " << this->get_ppm() << '\n';
    std::cout << "Smoothing length: " << this-> get_sm_len() << '\n';
    std::cout << "Particle mass: " << this->mass << '\n';
    std::cout << "Grid size: " << this->n_x << " x " << this->n_y << " x " << this->n_z << '\n';
    std::cout << "Number of blocks: " << n_blocks << '\n';
    std::cout << "Block size: " << this->size_blocks[0] << " x " << this->size_blocks[1] << " x " << this->size_blocks[2] << '\n';
}
