#include "grid.hpp"

Grid::Grid(int size_x, int size_y, int size_z) {
  this->size_x      = size_x;
  this->size_y      = size_y;
  this->size_z      = size_z;
  this->size        = size_x * size_y * size_z;
  this->n_particles = 0;
  this->blocks      = std::vector<Block>();
  for (int i = 0; i < this->size; ++i) { this->blocks.push_back(Block(i, size_x, size_y, size_z)); }
};

void Grid::populate(Simulation & sim, std::ifstream & inputFile) {
  int num_p = sim.get_num_p();
  for (int i = 0; i < num_p; i++) {
    Particle cur_particle = Particle(inputFile, i);
    int grid_index        = cur_particle.compute_grid_index(sim);
    this->blocks[grid_index].particles.push_back(cur_particle);
  }
};

void Grid::set_neighbors() {
  for (int block_id = 0; block_id < this->size; ++block_id) {
    Block & curr_block = this->blocks[block_id];
    // Iterate through neighbor blocks, taking care of edge cases using min and max
    // [!] We consider that a block is neighbor of itself (useful in computations)
    int x_start = std::max(curr_block.index_i - 1, 0);
    int x_end   = std::min(curr_block.index_i + 1, this->size_x - 1);
    int y_start = std::max(curr_block.index_j - 1, 0);
    int y_end   = std::min(curr_block.index_j + 1, this->size_y - 1);
    int z_start = std::max(curr_block.index_k - 1, 0);
    int z_end   = std::min(curr_block.index_k + 1, this->size_z - 1);
    for (int neig_id_x = x_start; neig_id_x <= x_end; ++neig_id_x) {
      for (int neig_id_y = y_start; neig_id_y <= y_end; ++neig_id_y) {
        for (int neig_id_z = z_start; neig_id_z <= z_end; ++neig_id_z) {
          // Compute global index of neighbor block
          int neig_id = neig_id_x + (this->size_x) * (neig_id_y + (this->size_y) * neig_id_z);
          // Get the neighbor block
          Block & neigh_block = this->blocks[neig_id];
          // Store pointers to neighbors in the current block
          curr_block.neighbours.push_back(&neigh_block);
        }
      }
    }
  }
}

bool Grid::cmp_trace(std::ifstream & trace) {
  // Read the total number of blocks from the binary file
  int32_t totalBlocks;
  trace.read(reinterpret_cast<char *>(&totalBlocks), sizeof(totalBlocks));
  // Check it matches size of the grid
  if (totalBlocks != static_cast<int32_t>(this->size)) {
    std::cerr << "Mismatch in the number of blocks in the grid." << std::endl;
    return false;
  }
  // Compare particle data block by block
  for (int i = 0; i < totalBlocks; ++i) {
    // Read the number of particles in the block from the binary file
    int64_t numParticles;
    trace.read(reinterpret_cast<char *>(&numParticles), sizeof(numParticles));

    // Check it matches the number of particles of the block in the trace
    if (numParticles != static_cast<int64_t>(this->blocks[i].particles.size())) {
      std::cerr << "Mismatch in the number of particles in block " << i << std::endl;
      std::cerr << "Grid: " << this->blocks[i].particles.size() << std::endl;
      std::cerr << "Expected: " << numParticles << std::endl;
      return false;
    }

    // Compare particle data field by field
    for (int j = 0; j < numParticles; ++j) {
      Particle particleFromFile;
      // trace.read(reinterpret_cast<char *>(&particleFromFile), sizeof(Particle));
      particleFromFile.pid     = read_binary_value<long>(trace);
      particleFromFile.posX    = read_binary_value<double>(trace);
      particleFromFile.posY    = read_binary_value<double>(trace);
      particleFromFile.posZ    = read_binary_value<double>(trace);
      particleFromFile.hvX     = read_binary_value<double>(trace);
      particleFromFile.hvY     = read_binary_value<double>(trace);
      particleFromFile.hvZ     = read_binary_value<double>(trace);
      particleFromFile.velX    = read_binary_value<double>(trace);
      particleFromFile.velY    = read_binary_value<double>(trace);
      particleFromFile.velZ    = read_binary_value<double>(trace);
      particleFromFile.density = read_binary_value<double>(trace);
      particleFromFile.accX    = read_binary_value<double>(trace);
      particleFromFile.accY    = read_binary_value<double>(trace);
      particleFromFile.accZ    = read_binary_value<double>(trace);

      Particle & particleFromVector = this->blocks[i].particles[j];

      // Compare each field of the particles
      double tolerance = 0;  // You can adjust the tolerance based on your specific use case

      if (fabs(particleFromFile.pid - particleFromVector.pid) > tolerance) {
        std::cout << "pid mismatch - Trace value: " << particleFromFile.pid
                  << ", Grid value: " << particleFromVector.pid << "\nBlock: " << i << std::endl;
        return false;
      }

      if (fabs(particleFromFile.posX - particleFromVector.posX) > tolerance) {
        std::cout << "posX mismatch - Trace value: " << particleFromFile.posX
                  << ", Grid value: " << particleFromVector.posX << "\nBlock: " << i << std::endl;
        std::cout << "Particle: " << j << std::endl;
        std::cout << "Difference: " << fabs(particleFromFile.posX - particleFromVector.posX)
                  << std::endl;
        return false;
      }

      if (fabs(particleFromFile.posY - particleFromVector.posY) > tolerance) {
        std::cout << "posY mismatch - Trace value: " << particleFromFile.posY
                  << ", Grid value: " << particleFromVector.posY << "\nBlock: " << i << std::endl;

        return false;
      }

      if (fabs(particleFromFile.posZ - particleFromVector.posZ) > tolerance) {
        std::cout << "posZ mismatch - Trace value: " << particleFromFile.posZ
                  << ", Grid value: " << particleFromVector.posZ << "\nBlock: " << i << std::endl;
        return false;
      }

      if (fabs(particleFromFile.hvX - particleFromVector.hvX) > tolerance) {
        std::cout << "hvX mismatch - Trace value: " << particleFromFile.hvX
                  << ", Grid value: " << particleFromVector.hvX << "\nBlock: " << i << std::endl;
        std::cout << "Difference: " << fabs(particleFromFile.hvX - particleFromVector.hvX)
                  << std::endl;
        return false;
      }

      if (fabs(particleFromFile.hvY - particleFromVector.hvY) > tolerance) {
        std::cout << "hvY mismatch - Trace value: " << particleFromFile.hvY
                  << ", Grid value: " << particleFromVector.hvY << "\nBlock: " << i << std::endl;
        return false;
      }

      if (fabs(particleFromFile.hvZ - particleFromVector.hvZ) > tolerance) {
        std::cout << "hvZ mismatch - Trace value: " << particleFromFile.hvZ
                  << ", Grid value: " << particleFromVector.hvZ << "\nBlock: " << i << std::endl;
        return false;
      }

      if (fabs(particleFromFile.velX - particleFromVector.velX) > tolerance) {
        std::cout << "velX mismatch - Trace value: " << particleFromFile.velX
                  << ", Grid value: " << particleFromVector.velX << "\nBlock: " << i << std::endl;
        std::cout << "Difference: " << fabs(particleFromFile.velX - particleFromVector.velX)
                  << std::endl;
        return false;
      }

      if (fabs(particleFromFile.velY - particleFromVector.velY) > tolerance) {
        std::cout << "velY mismatch - Trace value: " << particleFromFile.velY
                  << ", Grid value: " << particleFromVector.velY << "\nBlock: " << i << std::endl;
        return false;
      }

      if (fabs(particleFromFile.velZ - particleFromVector.velZ) > tolerance) {
        std::cout << "velZ mismatch - Trace value: " << particleFromFile.velZ
                  << ", Grid value: " << particleFromVector.velZ << "\nBlock: " << i << std::endl;
        std::cout << "Difference: " << fabs(particleFromFile.velZ - particleFromVector.velZ)
                  << std::endl;
        return false;
      }

      if (fabs(particleFromFile.density - particleFromVector.density) > tolerance) {
        std::cout << "density mismatch - Trace value: " << particleFromFile.density
                  << ", Grid value: " << particleFromVector.density << std::endl;
        std::cout << "Particle: " << j << " - Block: " << i << std::endl;
        std::cout << "Difference: " << fabs(particleFromFile.density - particleFromVector.density)
                  << "\n";
        return false;
      }

      if (fabs(particleFromFile.accX - particleFromVector.accX) > tolerance) {
        std::cout << "accX mismatch - Trace value: " << particleFromFile.accX
                  << ", Grid value: " << particleFromVector.accX << "\nBlock: " << i << std::endl;
        std::cout << "Particle: " << j << " - Block: " << i << std::endl;
        std::cout << "Difference: " << fabs(particleFromFile.accX - particleFromVector.accX)
                  << "\n";
        return false;
      }

      if (fabs(particleFromFile.accY - particleFromVector.accY) > tolerance) {
        std::cout << "accY mismatch - Trace value: " << particleFromFile.accY
                  << ", Grid value: " << particleFromVector.accY << "\n";
        std::cout << "Particle: " << j << " - Block: " << i << std::endl;
        std::cout << "Difference: " << fabs(particleFromFile.accY - particleFromVector.accY)
                  << std::endl;
        return false;
      }

      if (fabs(particleFromFile.accZ - particleFromVector.accZ) > tolerance) {
        std::cout << "accZ mismatch - Trace value: " << particleFromFile.accZ
                  << ", Grid value: " << particleFromVector.accZ << "\n";
        std::cout << "Particle: " << j << " - Block: " << i << std::endl;
        std::cout << "Difference: " << fabs(particleFromFile.accZ - particleFromVector.accZ)
                  << std::endl;
        return false;
      }
    }
  }

  std::cout << "Data comparison successful :)" << std::endl;
  return true;
}

// Only to fix precision issues
void Grid::set_to_trace(std::ifstream & trace) {
  // Read the total number of blocks from the binary file
  int32_t totalBlocks;
  trace.read(reinterpret_cast<char *>(&totalBlocks), sizeof(totalBlocks));

  // Compare particle data block by block
  for (int i = 0; i < totalBlocks; ++i) {
    // Read the number of particles in the block from the binary file
    int64_t numParticles;
    trace.read(reinterpret_cast<char *>(&numParticles), sizeof(numParticles));

    // Compare particle data field by field
    for (int j = 0; j < numParticles; ++j) {
      Particle particleFromFile;
      trace.read(reinterpret_cast<char *>(&particleFromFile), sizeof(Particle));

      Particle & particleFromVector = this->blocks[i].particles[j];

      if (particleFromFile.pid != particleFromVector.pid) {
        particleFromVector.pid = particleFromFile.pid;
      }

      if (particleFromFile.posX != particleFromVector.posX) {
        particleFromVector.posX = particleFromFile.posX;
      }

      if (particleFromFile.posY != particleFromVector.posY) {
        particleFromVector.posY = particleFromFile.posY;
      }

      if (particleFromFile.posZ != particleFromVector.posZ) {
        particleFromVector.posZ = particleFromFile.posZ;
      }

      if (particleFromFile.hvX != particleFromVector.hvX) {
        particleFromVector.hvX = particleFromFile.hvX;
      }

      if (particleFromFile.hvY != particleFromVector.hvY) {
        particleFromVector.hvY = particleFromFile.hvY;
      }

      if (particleFromFile.hvZ != particleFromVector.hvZ) {
        particleFromVector.hvZ = particleFromFile.hvZ;
      }

      if (particleFromFile.velX != particleFromVector.velX) {
        particleFromVector.velX = particleFromFile.velX;
      }

      if (particleFromFile.velY != particleFromVector.velY) {
        particleFromVector.velY = particleFromFile.velY;
      }

      if (particleFromFile.velZ != particleFromVector.velZ) {
        particleFromVector.velZ = particleFromFile.velZ;
      }

      if (particleFromFile.density != particleFromVector.density) {
        particleFromVector.density = particleFromFile.density;
      }

      if (particleFromFile.accX != particleFromVector.accX) {
        particleFromVector.accX = particleFromFile.accX;
      }

      if (particleFromFile.accY != particleFromVector.accY) {
        particleFromVector.accY = particleFromFile.accY;
      }

      if (particleFromFile.accZ != particleFromVector.accZ) {
        particleFromVector.accZ = particleFromFile.accZ;
      }
    }
  }

  std::cout << "Fixed grid!" << std::endl;
}

/*
void Grid::increase_all_dens(Simulation & sim) {
  std::set<std::pair<int, int>> proc_pairs;
  double const hSquared = std::pow(sim.get_sm_len(), 2);
  for (auto & block : this->blocks) {        // For each block in the grid
    for (auto & part_i : block.particles) {  // Iterate over particles within the current block
      for (auto & neigh_ptr : block.neighbours) {     // Iterate over neighbor block pointer
        for (auto & part_j : neigh_ptr->particles) {  // Iterate neigh particles
          int id_i = part_i.pid;
          int id_j = part_j.pid;
          if (id_i != id_j) {
            // (i,j) equivalent to (j,i)
            std::pair<int, int> particle_pair(std::min(id_i, id_j), std::max(id_i, id_j));
            // If not a processed pair
            if (proc_pairs.find(particle_pair) == proc_pairs.end()) {
              part_i.inc_part_dens(part_j, hSquared);
              proc_pairs.insert(particle_pair);  // Pair processed, keep track
            }
          }
        }
      }
    }
  }
} */

void Grid::increase_all_dens(Simulation & sim) {
  double const hSquared = std::pow(sim.get_sm_len(), 2);
  for (auto & block : this->blocks) {        // For each block in the grid
    for (auto & part_i : block.particles) {  // Iterate over particles within the current block
      for (auto & neigh_ptr : block.neighbours) {     // Iterate over neighbor block pointer
        for (auto & part_j : neigh_ptr->particles) {  // Iterate neigh particles
          int id_i = part_i.pid;
          int id_j = part_j.pid;
          if (id_j > id_i) {
            part_i.inc_part_dens(part_j, hSquared);
          }
        }
      }
    }
  }
}

void Grid::trans_all_dens(Simulation & sim) {
  double const sm_len         = sim.get_sm_len();
  double const h_pow6         = pow(sm_len, 6);
  double const prefactor_dens = (315 * sim.get_mass()) / (64 * std::numbers::pi * pow(sm_len, 9));
  // For each block in the grid
  for (auto & block : this->blocks) {
    // Iterate over particles within the current block
    for (auto & part_i : block.particles) {
      // Perform density transformation
      part_i.density = (part_i.density + h_pow6) * prefactor_dens;
    }
  }
}

// distanceSquared change
void Grid::increase_all_accs(Simulation & sim) {
  for (auto & block : this->blocks) {
    for (auto & part_i : block.particles) {
      for (auto & neigh_ptr : block.neighbours) {
        for (auto & part_j : neigh_ptr->particles) {
          int id_i = part_i.pid;
          int id_j = part_j.pid;
          if (id_j > id_i) {
            double distanceSquared = std::pow(part_i.posX - part_j.posX, 2) +
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

void Grid::part_collisions(Simulation & sim) {
  for (auto & block : this->blocks) {
    if (block.index_i == 0) {
      block.block_part_col_xmin(sim);
    } else if (block.index_i == this->size_x - 1) {
      block.block_part_col_xmax(sim);
    }
    if (block.index_j == 0) {
      block.block_part_col_ymin(sim);
    } else if (block.index_j == this->size_y - 1) {
      block.block_part_col_ymax(sim);
    }
    if (block.index_k == 0) {
      block.block_part_col_zmin(sim);
    } else if (block.index_k == this->size_z - 1) {
      block.block_part_col_zmax(sim);
    }
  }
}

void Grid::motion(Simulation & sim) {
  // For each block in the grid
  for (auto & block : this->blocks) {
    // For each particle in that block
    for (auto & particle : block.particles) {
      particle.posX = particle.posX + (particle.hvX * sim.delta_t) +
                      (particle.accX * sim.delta_t * sim.delta_t);
      particle.posY =
          particle.posY + particle.hvY * sim.delta_t + particle.accY * sim.delta_t * sim.delta_t;
      particle.posZ =
          particle.posZ + particle.hvZ * sim.delta_t + particle.accZ * sim.delta_t * sim.delta_t;

      particle.velX = particle.hvX + 0.5 * particle.accX * sim.delta_t;
      particle.velY = particle.hvY + 0.5 * particle.accY * sim.delta_t;
      particle.velZ = particle.hvZ + 0.5 * particle.accZ * sim.delta_t;

      particle.hvX = particle.hvX + particle.accX * sim.delta_t;
      particle.hvY = particle.hvY + particle.accY * sim.delta_t;
      particle.hvZ = particle.hvZ + particle.accZ * sim.delta_t;
    }
  }
}

void Grid::part_box_collisions(Simulation & sim) {
  for (auto & block : this->blocks) {
    if (block.index_i == 0) {
      block.boundint_xmin(sim);
    } else if (block.index_i == this->size_x - 1) {
      block.boundint_xmax(sim);
    }
    if (block.index_j == 0) {
      block.boundint_ymin(sim);
    } else if (block.index_j == this->size_y - 1) {
      block.boundint_ymax(sim);
    }
    if (block.index_k == 0) {
      block.boundint_zmin(sim);
    } else if (block.index_k == this->size_z - 1) {
      block.boundint_zmax(sim);
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
    int grid_index = par.compute_grid_index(sim);
    this->blocks[grid_index].particles.push_back(par);
  }
}

void Grid::init_acc() {
  for (auto & block : this->blocks) {
    for (auto & part : block.particles) {
      part.density = 0.0;
      part.accX    = 0.0;
      part.accY    = -9.8;
      part.accZ    = 0.0;
    }
  }
}

void Grid::gen_output(std::ofstream& out) {
  for (auto &block: this->blocks) {
    for (auto &particle: block.particles) {
      particle.write_particle_output(out);
    }
  }
}
