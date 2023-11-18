#include "block.hpp"

// Default constructor
Particle::Particle() { }

void Particle::inc_part_dens(Particle & part_j, double const hSquared) {
  double density_increase;
  // If not a processed pair
  double distanceSquared = std::pow(this->posX - part_j.posX, 2) +
                           std::pow(this->posY - part_j.posY, 2) +
                           std::pow(this->posZ - part_j.posZ, 2);
  if (distanceSquared < hSquared) {
    density_increase  = std::pow(hSquared - distanceSquared, 3);
    this->density    += density_increase;
    part_j.density   += density_increase;
  }
}

void Particle::inc_part_acc(Particle & part_j, Simulation & sim, double const distanceSquared) {
  double acc_inc_x, acc_inc_y, acc_inc_z;
  double const max_dis = 1e-12;
  double dist_ij       = sqrt(std::max(distanceSquared, max_dis));
  double fact_0_x      = this->posX - part_j.posX;
  double fact_0_y      = this->posY - part_j.posY;
  double fact_0_z      = this->posZ - part_j.posZ;
  double fact_4_x      = part_j.velX - this->velX;
  double fact_4_y      = part_j.velY - this->velY;
  double fact_4_z      = part_j.velZ - this->velZ;
  double fact_2        = (sim.get_sm_len() - dist_ij) * (sim.get_sm_len() - dist_ij) / dist_ij;
  double fact_3        = this->density + part_j.density - 2 * sim.fluid_density;
  acc_inc_x            = fact_0_x * sim.common_factor_acc * sim.fact_1_acc * fact_2 * fact_3 +
              fact_4_x * sim.common_factor_acc * sim.fact_5_acc;
  acc_inc_y = fact_0_y * sim.common_factor_acc * sim.fact_1_acc * fact_2 * fact_3 +
              fact_4_y * sim.common_factor_acc * sim.fact_5_acc;
  acc_inc_z = fact_0_z * sim.common_factor_acc * sim.fact_1_acc * fact_2 * fact_3 +
              fact_4_z * sim.common_factor_acc * sim.fact_5_acc;

  this->accX  += acc_inc_x / (this->density * part_j.density);
  this->accY  += acc_inc_y / (this->density * part_j.density);
  this->accZ  += acc_inc_z / (this->density * part_j.density);
  part_j.accX -= acc_inc_x / (this->density * part_j.density);
  part_j.accY -= acc_inc_y / (this->density * part_j.density);
  part_j.accZ -= acc_inc_z / (this->density * part_j.density);
}

// Particle constructor
Particle::Particle(std::ifstream & inputFile, int pid) {
  this->pid     = pid;
  this->posX    = read_binary_value<float>(inputFile);
  this->posY    = read_binary_value<float>(inputFile);
  this->posZ    = read_binary_value<float>(inputFile);
  this->hvX     = read_binary_value<float>(inputFile);
  this->hvY     = read_binary_value<float>(inputFile);
  this->hvZ     = read_binary_value<float>(inputFile);
  this->velX    = read_binary_value<float>(inputFile);
  this->velY    = read_binary_value<float>(inputFile);
  this->velZ    = read_binary_value<float>(inputFile);
  this->density = 0;
  this->accX    = 0;
  this->accY    = -9.8;
  this->accZ    = 0;
}

// Function to write Particle attributes in float precision
void Particle::write_particle_output(std::ofstream & outputFile) {
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

int Particle::compute_grid_index(Simulation & sim) {
  int i = floor((this->posX - sim.b_min[0]) / sim.size_blocks[0]);
  int j = floor((this->posY - sim.b_min[1]) / sim.size_blocks[1]);
  int k = floor((this->posZ - sim.b_min[2]) / sim.size_blocks[2]);

  if (i < 0) { i = 0; }
  if (i >= sim.n_x) { i = sim.n_x - 1; }
  if (j < 0) { j = 0; }
  if (j >= sim.n_y) { j = sim.n_y - 1; }
  if (k < 0) { k = 0; }
  if (k >= sim.n_z) { k = sim.n_z - 1; }

  int grid_index = i + sim.n_x * (j + sim.n_y * k);
  return grid_index;
}

Block::Block(int bid, int blocks_x, int blocks_y, int blocks_z) {
  this->bid        = bid;
  this->index_k    = bid / (blocks_x * blocks_y);
  int block_id_aux = bid % (blocks_x * blocks_y);
  this->index_j    = block_id_aux / blocks_x;
  this->index_i    = block_id_aux % blocks_x;
}

void Block::block_part_col_xmin(Simulation & sim) {
  double delta_coll;
  double const delta_coll_max = 1e-10;
  for (auto & particle : this->particles) {
    double temp_x = particle.posX + particle.hvX * sim.delta_t;
    delta_coll    = sim.d_p - (temp_x - sim.b_min[0]);
    if (delta_coll > delta_coll_max) {
      particle.accX += (sim.s_c * delta_coll - sim.d_v * particle.velX);
    }
  }
}

void Block::block_part_col_xmax(Simulation & sim) {
  double delta_coll;
  double const delta_coll_max = 1e-10;
  for (auto & particle : this->particles) {
    double temp_x = particle.posX + particle.hvX * sim.delta_t;
    delta_coll    = sim.d_p - (sim.b_max[0] - temp_x);
    if (delta_coll > delta_coll_max) {
      particle.accX -= (sim.s_c * delta_coll + sim.d_v * particle.velX);
    }
  }
}

void Block::block_part_col_ymin(Simulation & sim) {
  double delta_coll;
  double const delta_coll_max = 1e-10;
  for (auto & particle : this->particles) {
    double temp_y = particle.posY + particle.hvY * sim.delta_t;
    delta_coll    = sim.d_p - (temp_y - sim.b_min[1]);
    if (delta_coll > delta_coll_max) {
      particle.accY += (sim.s_c * delta_coll - sim.d_v * particle.velY);
    }
  }
}

void Block::block_part_col_ymax(Simulation & sim) {
  double delta_coll;
  double const delta_coll_max = 1e-10;
  for (auto & particle : this->particles) {
    double temp_y = particle.posY + particle.hvY * sim.delta_t;
    delta_coll    = sim.d_p - (sim.b_max[1] - temp_y);
    if (delta_coll > delta_coll_max) {
      particle.accY -= (sim.s_c * delta_coll + sim.d_v * particle.velY);
    }
  }
}

void Block::block_part_col_zmin(Simulation & sim) {
  double delta_coll;
  double const delta_coll_max = 1e-10;
  for (auto & particle : this->particles) {
    double temp_z = particle.posZ + particle.hvZ * sim.delta_t;
    delta_coll    = sim.d_p - (temp_z - sim.b_min[2]);
    if (delta_coll > delta_coll_max) {
      particle.accZ += (sim.s_c * delta_coll - sim.d_v * particle.velZ);
    }
  }
}

void Block::block_part_col_zmax(Simulation & sim) {
  double delta_coll;
  double const delta_coll_max = 1e-10;
  for (auto & particle : this->particles) {
    double temp_z = particle.posZ + particle.hvZ * sim.delta_t;
    delta_coll    = sim.d_p - (sim.b_max[2] - temp_z);
    if (delta_coll > delta_coll_max) {
      particle.accZ -= (sim.s_c * delta_coll + sim.d_v * particle.velZ);
    }
  }
}

void Block::boundint_xmin(Simulation & sim) {
  double dx;
  for (auto & particle : this->particles) {
    dx = particle.posX - sim.b_min[0];
    if (dx < 0) {
      particle.posX = sim.b_min[0] - dx;
      particle.velX = -particle.velX;
      particle.hvX  = -particle.hvX;
    }
  }
}

void Block::boundint_xmax(Simulation & sim) {
  double dx;
  for (auto & particle : this->particles) {
    dx = sim.b_max[0] - particle.posX;
    if (dx < 0) {
      particle.posX = sim.b_max[0] + dx;
      particle.velX = -particle.velX;
      particle.hvX  = -particle.hvX;
    }
  }
}

void Block::boundint_ymin(Simulation & sim) {
  double dy;
  for (auto & particle : this->particles) {
    dy = particle.posY - sim.b_min[1];
    if (dy < 0) {
      particle.posY = sim.b_min[1] - dy;
      particle.velY = -particle.velY;
      particle.hvY  = -particle.hvY;
    }
  }
}

void Block::boundint_ymax(Simulation & sim) {
  double dy;
  for (auto & particle : this->particles) {
    dy = sim.b_max[1] - particle.posY;
    if (dy < 0) {
      particle.posY = sim.b_max[1] + dy;
      particle.velY = -particle.velY;
      particle.hvY  = -particle.hvY;
    }
  }
}

void Block::boundint_zmin(Simulation & sim) {
  double dz;
  for (auto & particle : this->particles) {
    dz = particle.posZ - sim.b_min[2];
    if (dz < 0) {
      particle.posZ = sim.b_min[2] - dz;
      particle.velZ = -particle.velZ;
      particle.hvZ  = -particle.hvZ;
    }
  }
}

void Block::boundint_zmax(Simulation & sim) {
  double dz;
  for (auto & particle : this->particles) {
    dz = sim.b_max[2] - particle.posZ;
    if (dz < 0) {
      particle.posZ = sim.b_max[2] + dz;
      particle.velZ = -particle.velZ;
      particle.hvZ  = -particle.hvZ;
    }
  }
}
