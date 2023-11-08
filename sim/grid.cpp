#include "grid.hpp"

Grid::Grid(int size_x, int size_y, int size_z) {
    this->size_x = size_x;
    this->size_y = size_y;
    this->size_z = size_z;

    this->size = size_x * size_y * size_z;

    this->n_particles = 0;
    this->blocks = std::vector<Block>(this->size);
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
            if (particleFromFile.pid != particleFromVector.pid ||
                particleFromFile.posX != particleFromVector.posX ||
                particleFromFile.posY != particleFromVector.posY ||
                particleFromFile.posZ != particleFromVector.posZ ||
                // Compare other fields similarly...
                particleFromFile.accX != particleFromVector.accX ||
                particleFromFile.accY != particleFromVector.accY ||
                particleFromFile.accZ != particleFromVector.accZ) {
                std::cerr << "Mismatch in particle data in block " << i << ", particle " << j  << std::endl;
                return false;
            }
        }
    }

    std::cout << "Data comparison successful." << std::endl;
    return true;
}
