#include "grid.hpp"

Grid::Grid(int size_x, int size_y, int size_z) {
    this->size_x = size_x;
    this->size_y = size_y;
    this->size_z = size_z;

    this->size = size_x * size_y * size_z;

    this->n_particles = 0;

    this->blocks = std::vector<Block>();
    for (int i = 0; i < this->size; ++i) {
        this->blocks.push_back(Block(i));
    }
};

void Grid::populate(Simulation& sim, std::ifstream& inputFile) {
    int num_p = sim.get_num_p();
    for (int i = 0; i < num_p; i++) {
        Particle cur_particle = Particle(inputFile, i);
        int grid_index = cur_particle.compute_grid_index(sim);
        this->blocks[grid_index].particles.push_back(cur_particle);
    }
};

void Grid::display_grid(){
    int count = -1;
    for (auto block : this->blocks) {
        count++;
        std::cout << count << '\n';
        for (auto particle : block.particles) {
            std::cout << "ID: " << particle.pid  << " " << particle.posX << " " << particle.posY << " " << particle.posZ << std::endl;
        }
    }
}

void Grid::set_neighbors() {
    // For each block in the grid
    for (int block_id = 0; block_id < this->size; ++block_id) {
        // Retrieve indexes (i, j, k) from block_id
        int k        = block_id / (this->size_x * this->size_y);
        int block_id_aux = block_id % (this->size_x * this->size_y);
        int j        = block_id_aux / this->size_x;
        int i        = block_id_aux % this->size_x;
        // Get the current (actual) block
        Block& curr_block = this->blocks[block_id];
        // Iterate through neighbor blocks, taking care of edge cases using min (upper bounds) and max (lower bounds)
        // [!] We consider that a block is neighbor of itself (useful in computations)
        int x_start = std::max(i - 1, 0); 
        int x_end = std::min(i + 1, this->size_x - 1);
        int y_start = std::max(j - 1, 0);
        int y_end = std::min(j + 1, this->size_y - 1);
        int z_start = std::max(k - 1, 0);
        int z_end = std::min(k + 1, this->size_z - 1);
        for (int neig_id_x = x_start; neig_id_x <= x_end; ++neig_id_x) {
            for (int neig_id_y = y_start; neig_id_y <= y_end; ++neig_id_y) {
                for (int neig_id_z = z_start; neig_id_z <= z_end;++neig_id_z) {
                    // Compute global index of neighbor block
                    int neig_id = neig_id_x + (this->size_x) * (neig_id_y + (this->size_y) * neig_id_z);
                    // Get the neighbor block
                    Block& neigh_block = this->blocks[neig_id];
                    // Store pointers to neighbors in the current block
                    curr_block.neighbours.push_back(&neigh_block);
                }
            }
        }
    }
}

void Grid::test_neighbors()
{
    for (int i = 0; i < 5; ++i) {
        std::cout << "I am block " << i << ". These are my neighbors: " << std::endl;
        for (const auto neigh : this->blocks[i].neighbours) {
            std::cout << "Hey, I'm a neighbor. I am block: " << neigh->bid << std::endl;
        }
        std::cout << std::endl;
    }
}

bool Grid::cmp_trace(std::ifstream& trace)
{
    // Read the total number of blocks from the binary file
    int32_t totalBlocks;
    trace.read(reinterpret_cast<char*>(&totalBlocks), sizeof(totalBlocks));
    // Check it matches size of the grid
    if (totalBlocks != static_cast<int32_t>(this->size)) {
        std::cerr << "Mismatch in the number of blocks in the grid." << std::endl;
        return false;
    }
    // Compare particle data block by block
    for (int i = 0; i < totalBlocks; ++i) {
        // Read the number of particles in the block from the binary file
        int64_t numParticles;
        trace.read(reinterpret_cast<char*>(&numParticles), sizeof(numParticles));
        // Check it matches the number of particles of the block in the trace
        if (numParticles != static_cast<int64_t>(this->blocks[i].particles.size())) {
            std::cerr << "Mismatch in the number of particles in block " << i << std::endl;
            return false;
        }
        // Compare particle data field by field
        for (int j = 0; j < numParticles; ++j) {
            Particle particleFromFile;
            trace.read(reinterpret_cast<char*>(&particleFromFile), sizeof(Particle));

            Particle &particleFromVector = this->blocks[i].particles[j];

            // Compare each field of the particles
            double tolerance = 0; // You can adjust the tolerance based on your specific use case

            if (fabs(particleFromFile.pid - particleFromVector.pid) > tolerance) {
                std::cout << "pid mismatch - Trace value: " << particleFromFile.pid
                        << ", Grid value: " << particleFromVector.pid << "\nBlock: " << i << std::endl;
                return false;
            }

            if (fabs(particleFromFile.posX - particleFromVector.posX) > tolerance) {
                std::cout << "posX mismatch - Trace value: " << particleFromFile.posX
                        << ", Grid value: " << particleFromVector.posX << "\nBlock: " << i << std::endl;
                std::cout << "Particle: " << j << std::endl;
                std::cout << "Difference: " << fabs(particleFromFile.posX - particleFromVector.posX) << std::endl;
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
                std::cout << "Difference: " << fabs(particleFromFile.hvX - particleFromVector.hvX) << std::endl;
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
                std::cout << "Difference: " << fabs(particleFromFile.velX - particleFromVector.velX) << std::endl;
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
                std::cout << "Difference: " << fabs(particleFromFile.velZ - particleFromVector.velZ) << std::endl;
                return false;
            }

            if ((particleFromFile.density - particleFromVector.density) != 0) {
                std::cout << "density mismatch - Trace value: " << particleFromFile.density
                        << ", Grid value: " << particleFromVector.density << std::endl;
                std::cout << "Particle: " << j << " - Block: " << i << std::endl;
                std::cout << "Difference: " << fabs(particleFromFile.density - particleFromVector.density)  << "\n";
                this->blocks[i].particles[j].density = particleFromFile.density;
                return false;
            }

            if (fabs(particleFromFile.accX - particleFromVector.accX) > tolerance) {
                std::cout << "accX mismatch - Trace value: " << particleFromFile.accX
                        << ", Grid value: " << particleFromVector.accX << "\nBlock: " << i << std::endl;
                std::cout << "Particle: " << j << " - Block: " << i << std::endl;
                std::cout << "Difference: " << fabs(particleFromFile.accX - particleFromVector.accX)  << "\n";
                return false;
            }

            if (fabs(particleFromFile.accY - particleFromVector.accY) > tolerance) {
                std::cout << "accY mismatch - Trace value: " << particleFromFile.accY
                        << ", Grid value: " << particleFromVector.accY << "\nBlock: " << i << std::endl;
                return false;
            }

            if (fabs(particleFromFile.accZ - particleFromVector.accZ) > tolerance) {
                std::cout << "accZ mismatch - Trace value: " << particleFromFile.accZ
                        << ", Grid value: " << particleFromVector.accZ << "\nBlock: " << i << std::endl;
                return false;
            }
        }
    }

    std::cout << "Data comparison successful :)" << std::endl;
    return true;
}

// Only to fix precision issues
void Grid::set_to_trace(std::ifstream& trace){
    // Read the total number of blocks from the binary file
    int32_t totalBlocks;
    trace.read(reinterpret_cast<char*>(&totalBlocks), sizeof(totalBlocks));

    // Compare particle data block by block
    for (int i = 0; i < totalBlocks; ++i) {
        // Read the number of particles in the block from the binary file
        int64_t numParticles;
        trace.read(reinterpret_cast<char*>(&numParticles), sizeof(numParticles));

        // Compare particle data field by field
        for (int j = 0; j < numParticles; ++j) {
            Particle particleFromFile;
            trace.read(reinterpret_cast<char*>(&particleFromFile), sizeof(Particle));

            Particle &particleFromVector = this->blocks[i].particles[j];

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


// [!] TEMP: if successful, to be split in functions
void Grid::increase_all_dens(Simulation& sim) {

    std::set<std::pair<int, int>> proc_pairs;
    double sm_len = sim.get_sm_len();
    double const hSquared = sm_len * sm_len;
    double density_increase = 0;

    // For each block in the grid
    for (auto &block : this->blocks) {
        // Iterate over particles within the current block
        for (auto &part_i : block.particles) {
            // Iterate over neighbor block pointers
            for (auto &neigh_block_ptr : block.neighbours) {
                // Iterate over particles within the neighboring block
                for (auto& part_j : neigh_block_ptr->particles) {
                    int id_i = part_i.pid;
                    int id_j = part_j.pid;
                    // (i,j) equivalent to (j,i)
                    std::pair<int, int> particle_pair (std::min(id_i, id_j), std::max(id_i, id_j));
                    // If not a processed pair
                    if (proc_pairs.find(particle_pair) == proc_pairs.end() and id_i != id_j) {

                        // ------------------- DENSITY OPS -------------------
                        double diff_x = part_i.posX - part_j.posX;
                        double diff_y = part_i.posY - part_j.posY;
                        double diff_z = part_i.posZ - part_j.posZ;

                        double distanceSquared = (diff_x*diff_x) + (diff_y*diff_y) + (diff_z*diff_z);

                        if (distanceSquared < hSquared) {
                            double hMinusDist     = hSquared - distanceSquared;
                            density_increase      = hMinusDist * hMinusDist * hMinusDist;
                            part_i.density += density_increase;
                            part_j.density += density_increase;
                        }
                        // ------------------- DENSITY OPS -------------------

                        // Pair processed
                        proc_pairs.insert(particle_pair);
                    }
                }
            }
        }
    }
}

void Grid::trans_all_dens(Simulation &sim){

    double const sm_len = sim.get_sm_len();
    double const h_pow6 = pow(sm_len, 6);
    double const prefactor_dens = (315 * sim.get_mass()) / (64 * std::numbers::pi * pow(sm_len, 9));

    // For each block in the grid
    for (auto &block : this->blocks) { 
        // Iterate over particles within the current block
        for (auto &part_i : block.particles) {
             // Perform density transformation
            part_i.density = (part_i.density + h_pow6) * prefactor_dens;
        }
    }
}

void Grid::increase_all_accs(Simulation &sim) {
    std::set<std::pair<int, int>> proc_pairs;

    // Params and constants (future attributes)
    double sm_len = sim.get_sm_len();
    double const hSquared = sm_len * sm_len;
    double const f_dens = sim.fluid_density;
    double const common_factor = 45 / (std::numbers::pi*std::pow(sm_len,6));
    double const fact_1 = sim.get_mass()*sim.p_s / 2;
    double const fact_5 = sim.mew*sim.get_mass();

    // For each block in the grid
    for (auto &block : this->blocks) {
        // Iterate over particles within the current block
        for (auto &part_i : block.particles) {
            // Iterate over neighbor block pointers
            for (auto &neigh_block_ptr : block.neighbours) {
                // Iterate over particles within the neighboring block
                for (auto& part_j : neigh_block_ptr->particles) {

                    int id_i = part_i.pid;
                    int id_j = part_j.pid;
                    if (id_i != id_j) {   

                        double diff_x = part_i.posX - part_j.posX;
                        double diff_y = part_i.posY - part_j.posY;
                        double diff_z = part_i.posZ - part_j.posZ;
                        double distanceSquared = (diff_x*diff_x) + (diff_y*diff_y) + (diff_z*diff_z);

                        if (distanceSquared < hSquared) {
                            // Acc. computations
                            // (i,j) equivalent to (j,i)
                            std::pair<int, int> particle_pair (std::min(id_i, id_j), std::max(id_i, id_j));
                            // If not a processed pair
                            if (proc_pairs.find(particle_pair) == proc_pairs.end()) {
                                // ------------------- ACC OPS -------------------
                                double acc_inc_x = 0; 
                                double acc_inc_y = 0; 
                                double acc_inc_z = 0; 
                                
                                double const max_dis = 1e-12;
                                double dist_ij = sqrt(std::max(distanceSquared, max_dis));

                                double fact_0_x = part_i.posX - part_j.posX;
                                double fact_0_y = part_i.posY - part_j.posY;
                                double fact_0_z = part_i.posZ - part_j.posZ;

                                double fact_4_x = part_j.velX - part_i.velX;
                                double fact_4_y = part_j.velY - part_i.velY;
                                double fact_4_z = part_j.velZ - part_i.velZ;

                                double fact_2 = (sm_len-dist_ij)*(sm_len-dist_ij)/dist_ij;
                                double fact_3 = part_i.density+part_j.density - 2*f_dens;

                                acc_inc_x = fact_0_x*common_factor*fact_1*fact_2*fact_3 + fact_4_x*common_factor*fact_5;
                                acc_inc_y = fact_0_y*common_factor*fact_1*fact_2*fact_3 + fact_4_y*common_factor*fact_5;
                                acc_inc_z = fact_0_z*common_factor*fact_1*fact_2*fact_3 + fact_4_z*common_factor*fact_5;

                                // Update acc
                                part_i.accX += acc_inc_x / (part_i.density * part_j.density);
                                part_i.accY += acc_inc_y / (part_i.density * part_j.density);
                                part_i.accZ += acc_inc_z / (part_i.density * part_j.density);

                                part_j.accX -= acc_inc_x / (part_i.density * part_j.density);
                                part_j.accY -= acc_inc_y / (part_i.density * part_j.density);
                                part_j.accZ -= acc_inc_z / (part_i.density * part_j.density);
                                // ------------------- ACC OPS -------------------

                                // Pair processed
                                proc_pairs.insert(particle_pair);
                            }
                        }                            
                    }
                }
            }
        }
    }
}

void Grid::part_collisions(Simulation &sim) {

    double const delta_coll_max = 1e-10;
    double delta_coll;

    // For each block in the grid
    for (int block_id = 0; block_id < this->size; ++block_id) {
        // Retrieve indexes (i, j, k) from block_id
        int k        = block_id / (this->size_x * this->size_y);
        int block_id_aux = block_id % (this->size_x * this->size_y);
        int j        = block_id_aux / this->size_x;
        int i        = block_id_aux % this->size_x;

        // Boundary cases
        if (i == 0 || i == this->size_x - 1 || j == 0 || j == this->size_y - 1 || k == 0 || k == this->size_z - 1) {
        // For each particle in the current block
          for (auto & particle: this->blocks[block_id].particles) {
            if (i == 0) {
                
              particle.posX += particle.hvX * sim.delta_t;
              delta_coll    = sim.d_p - (particle.posX - sim.b_min[0]);
              if (delta_coll > delta_coll_max) {
                particle.accX += (sim.s_c * delta_coll - sim.d_v * particle.velX);
              }

            } else if (i == this->size_x - 1) {
              particle.posX += particle.hvX * sim.delta_t;
              delta_coll    = sim.d_p - (sim.b_max[0] - particle.posX);
              if (delta_coll > delta_coll_max) {
                particle.accX -= (sim.s_c * delta_coll + sim.d_v * particle.velX);
              }
            }

            if (j == 0) {
              particle.posY += particle.hvY * sim.delta_t;
              delta_coll    = sim.d_p - (particle.posY - sim.b_min[1]);
              if (delta_coll > delta_coll_max) {
                particle.accY += (sim.s_c * delta_coll - sim.d_v * particle.velY);
              }
            } else if (j == this->size_y - 1) {
              particle.posY += particle.hvY * sim.delta_t;
              delta_coll    = sim.d_p - (sim.b_max[1] - particle.posY);
              if (delta_coll > delta_coll_max) {
                particle.accY -= (sim.s_c * delta_coll + sim.d_v * particle.velY);
              }
            }

            if (k == 0) {
              particle.posZ += particle.hvZ * sim.delta_t;
              delta_coll    = sim.d_p - (particle.posZ - sim.b_min[2]);
              if (delta_coll > delta_coll_max) {
                particle.accZ += (sim.s_c * delta_coll - sim.d_v * particle.velZ);
              }
            } else if (k == this->size_z - 1) {
              particle.posZ += particle.hvZ * sim.delta_t;
              delta_coll    = sim.d_p - (sim.b_max[2] - particle.posZ);
              if (delta_coll > delta_coll_max) {
                particle.accZ -= (sim.s_c * delta_coll + sim.d_v * particle.velZ);
              }
            }
          }
        }
    }
}

void Grid::motion(Simulation &sim) {
    // For each block in the grid
    for (auto &block: this->blocks) {
        // For each particle in that block
        for (auto &particle: block.particles) {
            particle.posX = particle.posX + (particle.hvX * sim.delta_t) + (particle.accX * sim.delta_t * sim.delta_t);
            particle.posY = particle.posY + particle.hvY * sim.delta_t + particle.accY * sim.delta_t * sim.delta_t;
            particle.posZ = particle.posZ + particle.hvZ * sim.delta_t + particle.accZ * sim.delta_t * sim.delta_t;

            particle.velX = particle.hvX + 0.5 * particle.accX * sim.delta_t;
            particle.velY = particle.hvY + 0.5 * particle.accY * sim.delta_t;
            particle.velZ = particle.hvZ + 0.5 * particle.accZ * sim.delta_t;

            particle.hvX = particle.hvX + particle.accX * sim.delta_t;
            particle.hvY = particle.hvY + particle.accY * sim.delta_t;
            particle.hvZ = particle.hvZ + particle.accZ * sim.delta_t;
        }
    }
}

void Grid::part_box_collisions(Simulation &sim) {

    double const delta_coll_max = 1e-10;
    double delta_coll;
    double dx, dy, dz;

    // For each block in the grid
    for (int block_id = 0; block_id < this->size; ++block_id) {
        // Retrieve indexes (i, j, k) from block_id
        int k        = block_id / (this->size_x * this->size_y);
        int block_id_aux = block_id % (this->size_x * this->size_y);
        int j        = block_id_aux / this->size_x;
        int i        = block_id_aux % this->size_x;

        // Boundary cases
        if (i == 0 || i == this->size_x - 1 || j == 0 || j == this->size_y - 1 || k == 0 || k == this->size_z - 1) {
            // For each particle in the current block
            for (auto & particle: this->blocks[block_id].particles) {
                // Collisions with x-axis bounds
                if (i == 0) {
                    dx = particle.posX - sim.b_min[0];
                    if (dx < 0) {
                        particle.posX  = sim.b_min[0] - dx;
                        particle.velX  = -particle.velX;
                        particle.hvX = -particle.hvX;
                    }

                } if (i == this->size_x - 1) {

                    dx = sim.b_max[0] - particle.posX;
                    if (dx < 0) {
                        particle.posX  = sim.b_max[0] + dx;
                        particle.velX  = -particle.velX;
                        particle.hvX = -particle.hvX;
                    }
                }

                // Collisions with y-axis bounds
                if (j == 0) {
                    dy = particle.posY - sim.b_min[1];
                    if (dy < 0) {
                        particle.posY  = sim.b_min[1] - dy;
                        particle.velY  = -particle.velY;
                        particle.hvY = -particle.hvY;
                    }
                } if (j == this->size_y - 1) {
                    dy = sim.b_max[1] - particle.posY;
                    if (dy < 0) {
                        particle.posY  = sim.b_max[1] + dy;
                        particle.velY  = -particle.velY;
                        particle.hvY = -particle.hvY;
                    }
                }

                // Collisions with z-axis bounds
                if (k == 0) {
                    dz = particle.posZ - sim.b_min[2];
                    if (dz < 0) {
                        particle.posZ  = sim.b_min[2] - dz;
                        particle.velZ  = -particle.velZ;
                        particle.hvZ = -particle.hvZ;
                    }
                } if (k == this->size_z - 1) {
                    dz = sim.b_max[2] - particle.posZ;
                    if (dz < 0) {
                        particle.posZ  = sim.b_max[2] + dz;
                        particle.velZ  = -particle.velZ;
                        particle.hvZ = -particle.hvZ;
                    }
                }
            }
        }
    }
}

void Grid::repos(Simulation &sim) {
    // For each block in the grid
    for (int block_id = 0; block_id < this->size; ++block_id) {
        // For each particle in that block
        for (auto & curr_par : this->blocks[block_id].particles) {
            int grid_index = curr_par.compute_grid_index(sim);
            if (block_id != grid_index) {
                Particle cp_par = curr_par;
                this->blocks[grid_index].particles.push_back(cp_par);
                curr_par.pid = -1; // Keep track of pars to delete from block block_id
            }
        }
        // Erase particles with pid = -1 from this block
        auto remove_end = std::remove_if(this->blocks[block_id].particles.begin(), this->blocks[block_id].particles.end(), [](const Particle& particle) {
            return particle.pid == -1;
        });
        this->blocks[block_id].particles.erase(remove_end, this->blocks[block_id].particles.end());
    }
}