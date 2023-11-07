#include "block.hpp"

// Particle constructor
Particle::Particle(std::ifstream& inputFile, int pid){
    this->pid = pid;
    posX = read_binary_value<float>(inputFile);
    posY = read_binary_value<float>(inputFile);
    posZ = read_binary_value<float>(inputFile);
    hvX = read_binary_value<float>(inputFile);
    hvY = read_binary_value<float>(inputFile);
    hvZ = read_binary_value<float>(inputFile);
    velX = read_binary_value<float>(inputFile);
    velY = read_binary_value<float>(inputFile);
    velZ = read_binary_value<float>(inputFile);
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