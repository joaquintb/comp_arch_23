#include "block.hpp"

// Default constructor
Particle::Particle() = default;

void Particle::inc_part_dens(Particle & part_j, double const hSquared) {
  double const distanceSquared = std::pow(this->posX - part_j.posX, 2) +
                                 std::pow(this->posY - part_j.posY, 2) +
                                 std::pow(this->posZ - part_j.posZ, 2);
  if (distanceSquared < hSquared) {
    double const density_increase  = std::pow(hSquared - distanceSquared, 3);
    this->density                 += density_increase;
    part_j.density                += density_increase;
  }
}

void Particle::inc_part_acc(Particle & part_j, Simulation & sim, double const distanceSquared) {
  double const max_dis   = 1e-12;
  double const dist_ij   = sqrt(std::max(distanceSquared, max_dis));
  double const fact_0_x  = this->posX - part_j.posX;
  double const fact_0_y  = this->posY - part_j.posY;
  double const fact_0_z  = this->posZ - part_j.posZ;
  double const fact_4_x  = part_j.velX - this->velX;
  double const fact_4_y  = part_j.velY - this->velY;
  double const fact_4_z  = part_j.velZ - this->velZ;
  double const fact_2    = (sim.sm_len - dist_ij) * (sim.sm_len - dist_ij) / dist_ij;
  double const fact_3    = this->density + part_j.density - 2 * Simulation::fluid_density;
  double const acc_inc_x = fact_0_x * sim.common_factor_acc * sim.fact_1_acc * fact_2 * fact_3 +
                           fact_4_x * sim.common_factor_acc * sim.fact_5_acc;
  double const acc_inc_y = fact_0_y * sim.common_factor_acc * sim.fact_1_acc * fact_2 * fact_3 +
                           fact_4_y * sim.common_factor_acc * sim.fact_5_acc;
  double const acc_inc_z = fact_0_z * sim.common_factor_acc * sim.fact_1_acc * fact_2 * fact_3 +
                           fact_4_z * sim.common_factor_acc * sim.fact_5_acc;
  this->accX  += acc_inc_x / (this->density * part_j.density);
  this->accY  += acc_inc_y / (this->density * part_j.density);
  this->accZ  += acc_inc_z / (this->density * part_j.density);
  part_j.accX -= acc_inc_x / (this->density * part_j.density);
  part_j.accY -= acc_inc_y / (this->density * part_j.density);
  part_j.accZ -= acc_inc_z / (this->density * part_j.density);
}

// Particle constructor
Particle::Particle(std::ifstream & inputFile, int pid)
  : pid(pid), posX(read_binary_value<float>(inputFile)), posY(read_binary_value<float>(inputFile)),
    posZ(read_binary_value<float>(inputFile)), hvX(read_binary_value<float>(inputFile)),
    hvY(read_binary_value<float>(inputFile)), hvZ(read_binary_value<float>(inputFile)),
    velX(read_binary_value<float>(inputFile)), velY(read_binary_value<float>(inputFile)),
    velZ(read_binary_value<float>(inputFile)), accY(Simulation::gravity) { }

// Function to write Particle attributes in float precision
void Particle::write_particle_output(std::ofstream & outputFile) const {
  write_binary_value(static_cast<float>(this->posX), outputFile);
  write_binary_value(static_cast<float>(this->posY), outputFile);
  write_binary_value(static_cast<float>(this->posZ), outputFile);
  write_binary_value(static_cast<float>(this->hvX), outputFile);
  write_binary_value(static_cast<float>(this->hvY), outputFile);
  write_binary_value(static_cast<float>(this->hvZ), outputFile);
  write_binary_value(static_cast<float>(this->velX), outputFile);
  write_binary_value(static_cast<float>(this->velY), outputFile);
  write_binary_value(static_cast<float>(this->velZ), outputFile);
}

int Particle::compute_grid_index(Simulation & sim) const {
  int index_i = floor((this->posX - Simulation::b_min[0]) / sim.size_blocks[0]);
  int index_j = floor((this->posY - Simulation::b_min[1]) / sim.size_blocks[1]);
  int index_k = floor((this->posZ - Simulation::b_min[2]) / sim.size_blocks[2]);

  if (index_i < 0) { index_i = 0; }
  if (index_i >= sim.n_x) { index_i = sim.n_x - 1; }
  if (index_j < 0) { index_j = 0; }
  if (index_j >= sim.n_y) { index_j = sim.n_y - 1; }
  if (index_k < 0) { index_k = 0; }
  if (index_k >= sim.n_z) { index_k = sim.n_z - 1; }

  int const grid_index = index_i + sim.n_x * (index_j + sim.n_y * index_k);
  return grid_index;
}

Block::Block(int bid, int blocks_x, int blocks_y) : bid(bid), index_k(bid / (blocks_x * blocks_y)) {
  int const block_id_aux = bid % (blocks_x * blocks_y);
  this->index_j          = block_id_aux / blocks_x;
  this->index_i          = block_id_aux % blocks_x;
}

void Block::block_part_col_xmin() {
  for (auto & particle : this->particles) {
    double const temp_x     = particle.posX + particle.hvX * Simulation::delta_t;
    double const delta_coll = Simulation::d_p - (temp_x - Simulation::b_min[0]);
    if (delta_coll > Simulation::delta_coll_max) {
      particle.accX += (Simulation::s_c * delta_coll - Simulation::d_v * particle.velX);
    }
  }
}

void Block::block_part_col_xmax() {
  for (auto & particle : this->particles) {
    double const temp_x     = particle.posX + particle.hvX * Simulation::delta_t;
    double const delta_coll = Simulation::d_p - (Simulation::b_max[0] - temp_x);
    if (delta_coll > Simulation::delta_coll_max) {
      particle.accX -= (Simulation::s_c * delta_coll + Simulation::d_v * particle.velX);
    }
  }
}

void Block::block_part_col_ymin() {
  for (auto & particle : this->particles) {
    double const temp_y     = particle.posY + particle.hvY * Simulation::delta_t;
    double const delta_coll = Simulation::d_p - (temp_y - Simulation::b_min[1]);
    if (delta_coll > Simulation::delta_coll_max) {
      particle.accY += (Simulation::s_c * delta_coll - Simulation::d_v * particle.velY);
    }
  }
}

void Block::block_part_col_ymax() {
  for (auto & particle : this->particles) {
    double const temp_y     = particle.posY + particle.hvY * Simulation::delta_t;
    double const delta_coll = Simulation::d_p - (Simulation::b_max[1] - temp_y);
    if (delta_coll > Simulation::delta_coll_max) {
      particle.accY -= (Simulation::s_c * delta_coll + Simulation::d_v * particle.velY);
    }
  }
}

void Block::block_part_col_zmin() {
  for (auto & particle : this->particles) {
    double const temp_z     = particle.posZ + particle.hvZ * Simulation::delta_t;
    double const delta_coll = Simulation::d_p - (temp_z - Simulation::b_min[2]);
    if (delta_coll > Simulation::delta_coll_max) {
      particle.accZ += (Simulation::s_c * delta_coll - Simulation::d_v * particle.velZ);
    }
  }
}

void Block::block_part_col_zmax() {
  for (auto & particle : this->particles) {
    double const temp_z     = particle.posZ + particle.hvZ * Simulation::delta_t;
    double const delta_coll = Simulation::d_p - (Simulation::b_max[2] - temp_z);
    if (delta_coll > Simulation::delta_coll_max) {
      particle.accZ -= (Simulation::s_c * delta_coll + Simulation::d_v * particle.velZ);
    }
  }
}

void Block::boundint_xmin() {
  for (auto & particle : this->particles) {
    double const d_x = particle.posX - Simulation::b_min[0];
    if (d_x < 0) {
      particle.posX = Simulation::b_min[0] - d_x;
      particle.velX = -particle.velX;
      particle.hvX  = -particle.hvX;
    }
  }
}

void Block::boundint_xmax() {
  for (auto & particle : this->particles) {
    double const d_x = Simulation::b_max[0] - particle.posX;
    if (d_x < 0) {
      particle.posX = Simulation::b_max[0] + d_x;
      particle.velX = -particle.velX;
      particle.hvX  = -particle.hvX;
    }
  }
}

void Block::boundint_ymin() {
  for (auto & particle : this->particles) {
    double const d_y = particle.posY - Simulation::b_min[1];
    if (d_y < 0) {
      particle.posY = Simulation::b_min[1] - d_y;
      particle.velY = -particle.velY;
      particle.hvY  = -particle.hvY;
    }
  }
}

void Block::boundint_ymax() {
  for (auto & particle : this->particles) {
    double const d_y = Simulation::b_max[1] - particle.posY;
    if (d_y < 0) {
      particle.posY = Simulation::b_max[1] + d_y;
      particle.velY = -particle.velY;
      particle.hvY  = -particle.hvY;
    }
  }
}

void Block::boundint_zmin() {
  for (auto & particle : this->particles) {
    double const d_z = particle.posZ - Simulation::b_min[2];
    if (d_z < 0) {
      particle.posZ = Simulation::b_min[2] - d_z;
      particle.velZ = -particle.velZ;
      particle.hvZ  = -particle.hvZ;
    }
  }
}

void Block::boundint_zmax() {
  for (auto & particle : this->particles) {
    double const d_z = Simulation::b_max[2] - particle.posZ;
    if (d_z < 0) {
      particle.posZ = Simulation::b_max[2] + d_z;
      particle.velZ = -particle.velZ;
      particle.hvZ  = -particle.hvZ;
    }
  }
}
