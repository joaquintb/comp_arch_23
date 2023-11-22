#include "grid.hpp"

Grid::Grid(int size_x, int size_y, int size_z)
  : size_x(size_x), size_y(size_y), size_z(size_z), size(size_x * size_y * size_z) {
  this->blocks = std::vector<Block>();
  for (int i = 0; i < this->size; ++i) { this->blocks.emplace_back(i, size_x, size_y); }
}

void Grid::populate(Simulation & sim, std::ifstream & inputFile) {
  int num_p = 0;
  while (inputFile.peek() != EOF) {
    Particle const cur_particle = Particle(inputFile, num_p);
    int const grid_index        = cur_particle.compute_grid_index(sim);
    this->blocks[grid_index].particles.push_back(cur_particle);
    num_p++;
  }
  if (num_p != sim.num_p) {
    std::cerr << "Error: Number of particles mismatch. Header: " << sim.num_p
              << ", Found: " << num_p << ".\n";
    std::exit(Simulation::particle_error_code);
  }
}

void Grid::set_neighbors() {
  for (int block_id = 0; block_id < this->size; ++block_id) {
    Block & curr_block = this->blocks[block_id];
    // Iterate through neighbor blocks, taking care of edge cases using min and max
    // [!] We consider that a block is neighbor of itself (useful in computations)
    int const x_start = std::max(curr_block.index_i - 1, 0);
    int const x_end   = std::min(curr_block.index_i + 1, this->size_x - 1);
    int const y_start = std::max(curr_block.index_j - 1, 0);
    int const y_end   = std::min(curr_block.index_j + 1, this->size_y - 1);
    int const z_start = std::max(curr_block.index_k - 1, 0);
    int const z_end   = std::min(curr_block.index_k + 1, this->size_z - 1);
    for (int neig_id_x = x_start; neig_id_x <= x_end; ++neig_id_x) {
      for (int neig_id_y = y_start; neig_id_y <= y_end; ++neig_id_y) {
        for (int neig_id_z = z_start; neig_id_z <= z_end; ++neig_id_z) {
          // Compute global index of neighbor block
          int const neig_id = neig_id_x + (this->size_x) * (neig_id_y + (this->size_y) * neig_id_z);
          // Get the neighbor block
          Block & neigh_block = this->blocks[neig_id];
          // Store pointers to neighbors in the current block
          curr_block.neighbours.push_back(&neigh_block);
        }
      }
    }
  }
}

void Grid::increase_all_dens(Simulation & sim) {
  double const hSquared = sim.sm_len_sq;
  for (auto & block : this->blocks) {        // For each block in the grid
    for (auto & part_i : block.particles) {  // Iterate over particles within the current block
      for (auto & neigh_ptr : block.neighbours) {     // Iterate over neighbor block pointer
        for (auto & part_j : neigh_ptr->particles) {  // Iterate neigh particles
          long const id_i = part_i.pid;
          long const id_j = part_j.pid;
          if (id_j > id_i) { part_i.inc_part_dens(part_j, hSquared); }
        }
      }
    }
  }
}

void Grid::trans_all_dens(Simulation & sim) {
  // For each block in the grid
  for (auto & block : this->blocks) {
    // Iterate over particles within the current block
    for (auto & part_i : block.particles) {
      // Perform density transformation
      part_i.density = (part_i.density + sim.sm_len_six) * sim.prefactor_dens;
    }
  }
}

// distanceSquared change
void Grid::increase_all_accs(Simulation & sim) {
  for (auto & block : this->blocks) {
    for (auto & part_i : block.particles) {
      for (auto & neigh_ptr : block.neighbours) {
        for (auto & part_j : neigh_ptr->particles) {
          long const id_i = part_i.pid;
          long const id_j = part_j.pid;
          if (id_j > id_i) {
            double const distanceSquared = std::pow(part_i.posX - part_j.posX, 2) +
                                           std::pow(part_i.posY - part_j.posY, 2) +
                                           std::pow(part_i.posZ - part_j.posZ, 2);
            if (distanceSquared < sim.sm_len_sq) {
              part_i.inc_part_acc(part_j, sim, distanceSquared);
            }
          }
        }
      }
    }
  }
}

void Grid::part_collisions() {
  for (auto & block : this->blocks) {
    if (block.index_i == 0) {
      block.block_part_col_xmin();
    } else if (block.index_i == this->size_x - 1) {
      block.block_part_col_xmax();
    }
    if (block.index_j == 0) {
      block.block_part_col_ymin();
    } else if (block.index_j == this->size_y - 1) {
      block.block_part_col_ymax();
    }
    if (block.index_k == 0) {
      block.block_part_col_zmin();
    } else if (block.index_k == this->size_z - 1) {
      block.block_part_col_zmax();
    }
  }
}

void Grid::motion() {
  double const half = 0.5;
  // For each block in the grid
  for (auto & block : this->blocks) {
    // For each particle in that block
    for (auto & particle : block.particles) {
      particle.posX = particle.posX + particle.hvX * Simulation::delta_t +
                      particle.accX * Simulation::delta_t * Simulation::delta_t;
      particle.posY = particle.posY + particle.hvY * Simulation::delta_t +
                      particle.accY * Simulation::delta_t * Simulation::delta_t;
      particle.posZ = particle.posZ + particle.hvZ * Simulation::delta_t +
                      particle.accZ * Simulation::delta_t * Simulation::delta_t;

      particle.velX = particle.hvX + half * particle.accX * Simulation::delta_t;
      particle.velY = particle.hvY + half * particle.accY * Simulation::delta_t;
      particle.velZ = particle.hvZ + half * particle.accZ * Simulation::delta_t;

      particle.hvX = particle.hvX + particle.accX * Simulation::delta_t;
      particle.hvY = particle.hvY + particle.accY * Simulation::delta_t;
      particle.hvZ = particle.hvZ + particle.accZ * Simulation::delta_t;
    }
  }
}

void Grid::part_box_collisions() {
  for (auto & block : this->blocks) {
    if (block.index_i == 0) {
      block.boundint_xmin();
    } else if (block.index_i == this->size_x - 1) {
      block.boundint_xmax();
    }
    if (block.index_j == 0) {
      block.boundint_ymin();
    } else if (block.index_j == this->size_y - 1) {
      block.boundint_ymax();
    }
    if (block.index_k == 0) {
      block.boundint_zmin();
    } else if (block.index_k == this->size_z - 1) {
      block.boundint_zmax();
    }
  }
}

void Grid::repos(Simulation & sim) {
  std::vector<Particle> all_part;

  for (auto & block : this->blocks) {
    for (auto & par : block.particles) { all_part.push_back(par); }
  }

  for (auto & block : this->blocks) { block.particles.clear(); }

  for (auto & par : all_part) {
    int const grid_index = par.compute_grid_index(sim);
    this->blocks[grid_index].particles.push_back(par);
  }
}

void Grid::init_acc() {
  for (auto & block : this->blocks) {
    for (auto & part : block.particles) {
      part.density = 0.0;
      part.accX    = 0.0;
      part.accY    = Simulation::gravity;
      part.accZ    = 0.0;
    }
  }
}

void Grid::gen_output(std::ofstream & out) {
  // Collect all particles from all blocks into a single vector
  std::vector<Particle> all_particles;
  for (auto & block : this->blocks) {
    all_particles.insert(all_particles.end(), block.particles.begin(), block.particles.end());
  }

  // Sort all particles based on pid attribute
  std::sort(all_particles.begin(), all_particles.end(),
            [](Particle const & part_1, Particle const & part_2) {
              return part_1.pid < part_2.pid;
            });

  // Write the sorted particles to the output
  for (auto & particle : all_particles) { particle.write_particle_output(out); }
}

void Grid::simulate(int n_steps, Simulation & sim) {
  for (int time_step = 0; time_step < n_steps; ++time_step) {
    if (time_step > 0) {
      this->repos(sim);
      this->init_acc();
    }

    this->increase_all_dens(sim);
    this->trans_all_dens(sim);
    this->increase_all_accs(sim);
    this->part_collisions();
    this->motion();
    this->part_box_collisions();
  }
}
