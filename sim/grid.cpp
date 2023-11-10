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
    for (auto block : this->blocks) {
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

            const Particle& particleFromVector = this->blocks[i].particles[j];

            // Compare each field of the particles
            double tolerance = 0.000000000001; // You can adjust the tolerance based on your specific use case

            if (fabs(particleFromFile.density - particleFromVector.density) > tolerance)  {
                std::cout << "Trace value: " << particleFromFile.density << std::endl;
                std::cout << "Grid value: " << particleFromVector.density << std::endl;
                std::cerr << "Mismatch in particle data in block " << i << ", particle " << j << std::endl;
                return false;
            }
        }
    }

    std::cout << "Data comparison successful." << std::endl;
    return true;
}


// [!] TEMP: if successful, to be split in functions
void Grid::increase_all_dens(Simulation& sim) {

    std::set<std::pair<int, int>> proc_pairs;

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
                        double density_increase = 0;
                        double diff_x = part_i.posX - part_j.posX;
                        double diff_y = part_i.posY - part_j.posY;
                        double diff_z = part_i.posZ - part_j.posZ;

                        double distanceSquared = (diff_x)*(diff_x) + (diff_y)*(diff_y) + (diff_z)*(diff_z);

                        double h = sim.get_sm_len();
                        double const hSquared = h * h;

                        if (distanceSquared < hSquared) {
                            double hMinusDist     = hSquared - distanceSquared;
                            density_increase      = hMinusDist * hMinusDist * hMinusDist;
                        }

                        part_i.density += density_increase;
                        part_j.density += density_increase;

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
