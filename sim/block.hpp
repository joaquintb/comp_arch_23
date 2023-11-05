//
// Created by jtorresb on 5/11/23.
//
#include <vector>
#include <fstream>
// #include "read_binary.hpp"

#ifndef COMP_ARCH_23_BLOCK_HPP
#define COMP_ARCH_23_BLOCK_HPP


struct Particle{
    double px, py, pz, hvx, hvy, hvz, vx, vy, vz, density, ax, ay, az;
    int id;

    Particle(std::ifstream& file, int id);
    template <typename T>
    T read_binary_value(std::istream &);
    template <typename T>
    void write_binary_value(T, std::ostream & ); 
    template <typename T>
    char const * as_buffer(T const &);
    template <typename T>
    char * as_writable_buffer(T &);
};

class Block{
    // Attributes 
    std::vector<Particle> particles;
    std::vector<Block*> neighbours; // max size of 26

    // Methods 
    void reposition_particles();
    void compute_forces(); 
    void process_collisions();
    void move_particles(); 
    void process_boundaries();
};
#endif //COMP_ARCH_23_BLOCK_HPP
