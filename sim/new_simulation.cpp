#include "simulation.hpp"

Simulation::Simulation(float ppm, int num_p)
  : ppm(static_cast<double>(ppm)), sm_len(Simulation::radius / this->ppm),
    mass(Simulation::fluid_density / std::pow(this->ppm, 3)), num_p(num_p),
    n_x(floor((Simulation::b_max[0] - Simulation::b_min[0]) / this->sm_len)),
    n_y(floor((Simulation::b_max[1] - Simulation::b_min[1]) / this->sm_len)),
    n_z(floor((Simulation::b_max[2] - Simulation::b_min[2]) / this->sm_len)),
    num_blocks(this->n_x * this->n_y * this->n_z), sm_len_sq(std::pow(this->sm_len, 2)),
    fact_5_acc(Simulation::mew * this->mass) {

  this->size_blocks       = {(Simulation::b_max[0] - Simulation::b_min[0]) / this->n_x,
                             (Simulation::b_max[1] - Simulation::b_min[1]) / this->n_y,
                             (Simulation::b_max[2] - Simulation::b_min[2]) / this->n_z};
  double const exp        = 6.0;
  double const half       = 0.5;
  this->common_factor_acc = forty_five_over_pi / std::pow(this->sm_len, exp);
  this->fact_1_acc        = this->mass * Simulation::p_s * half;
}

void Simulation::check_positive_particles() const {
  const int err_num_p = -5;
  if(this->num_p <= 0){
    std::cerr << "Error: Invalid number of particles: " << this->num_p << ".\n";
    std::exit(err_num_p);
  }
}

void Simulation::print_sim_values() {
  std::cout << "Number of particles: " << this->num_p<< '\n';
  std::cout << "Particles per meter: " << this->ppm << '\n';
  std::cout << "Smoothing length: " << this->sm_len << '\n';
  std::cout << "Particle mass: " << this->mass << '\n';
  std::cout << "Grid size: " << this->n_x << " x " << this->n_y << " x " << this->n_z << '\n';
  std::cout << "Number of blocks: " << this->num_blocks << '\n';
  std::cout << "Block size: " << this->size_blocks[0] << " x " << this->size_blocks[1] << " x " <<
      this->size_blocks[2] << '\n';
}