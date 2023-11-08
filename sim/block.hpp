#ifndef COMP_ARCH_23_BLOCK_HPP
#define COMP_ARCH_23_BLOCK_HPP

#include <vector>
#include <fstream>

#include "simulation.hpp"

struct Particle{
    int64_t pid;
    double posX, posY, posZ;
    double hvX, hvY, hvZ;
    double velX, velY, velZ;
    double density;
    double accX, accY, accZ;

    Particle();

    Particle(std::ifstream& inputFile, int pid);

    int compute_grid_index(Simulation& sim);

    void write_particle_trace(std::ofstream& outputFile);

    template <typename T>
    requires(std::is_integral_v<T> or std::is_floating_point_v<T>)
    T read_binary_value(std::istream & is) {
        T value{};
        is.read(as_writable_buffer(value), sizeof(value));
        return value;
    }

    template <typename T>
    requires(std::is_integral_v<T> or std::is_floating_point_v<T>)
    void write_binary_value(T value, std::ostream & os) {
        os.write(as_buffer(value), sizeof(value));
    }

    template <typename T>
    requires(std::is_integral_v<T> or std::is_floating_point_v<T>)
    char const * as_buffer(T const & value) {
        // NOLINTNEXTLINE(cppcoreguidelines-pro-type-reinterpret-cast)
        return reinterpret_cast<char const *>(&value);
    }

    template <typename T>
    char * as_writable_buffer(T & value) {
        // NOLINTNEXTLINE(cppcoreguidelines-pro-type-reinterpret-cast)
        return reinterpret_cast<char *>(&value);
    }
};

class Block{
    // Attributes 
public:
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
