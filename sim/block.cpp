#include "block.hpp"

// Default constructor
Particle::Particle(){
}

void Particle::inc_part_dens(Particle & part_j, double const hSquared) {
    double density_increase;
    // If not a processed pair
    double distanceSquared = std::pow(this->posX - part_j.posX, 2) +
                             std::pow(this->posY - part_j.posY, 2) +
                             std::pow(this->posZ - part_j.posZ, 2);
    if (distanceSquared < hSquared) {
      density_increase  = std::pow(hSquared - distanceSquared, 3);
      this->density   += density_increase;
      part_j.density   += density_increase;
    }
}

// Particle constructor
Particle::Particle(std::ifstream& inputFile, int pid) {
    this->pid = pid;
    this->posX = read_binary_value<float>(inputFile);
    this->posY = read_binary_value<float>(inputFile);
    this->posZ = read_binary_value<float>(inputFile);
    this->hvX = read_binary_value<float>(inputFile);
    this->hvY = read_binary_value<float>(inputFile);
    this->hvZ = read_binary_value<float>(inputFile);
    this->velX = read_binary_value<float>(inputFile);
    this->velY = read_binary_value<float>(inputFile);
    this->velZ = read_binary_value<float>(inputFile);
    this->density = 0;
    this->accX = 0;
    this->accY = -9.8;
    this->accZ = 0;
}

// Function to write Particle attributes in double precision
void Particle::write_particle_trace(std::ofstream& outputFile) {
    this->write_binary_value(this->pid, outputFile);
    this->write_binary_value(this->posX, outputFile);
    this->write_binary_value(this->posY, outputFile);
    this->write_binary_value(this->posZ, outputFile);
    this->write_binary_value(this->hvX, outputFile);
    this->write_binary_value(this->hvY, outputFile);
    this->write_binary_value(this->hvZ, outputFile);
    this->write_binary_value(this->velX, outputFile);
    this->write_binary_value(this->velY, outputFile);
    this->write_binary_value(this->velZ, outputFile);
    this->write_binary_value(this->density, outputFile);
    this->write_binary_value(this->accX, outputFile);
    this->write_binary_value(this->accY, outputFile);
    this->write_binary_value(this->accZ, outputFile);
}

int Particle::compute_grid_index(Simulation& sim) {
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

Block::Block (int bid, int blocks_x, int blocks_y, int blocks_z){
    this->bid = bid;
    this->index_k  = bid / (blocks_x * blocks_y);
    int block_id_aux = bid % (blocks_x * blocks_y);
    this->index_j  = block_id_aux / blocks_x;
    this->index_i    = block_id_aux % blocks_x;
}
