#include "block.hpp"

// Particle constructor
Particle::Particle(std::ifstream& inputFile, int pid){
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