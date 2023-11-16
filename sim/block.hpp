#ifndef COMP_ARCH_23_BLOCK_HPP
#define COMP_ARCH_23_BLOCK_HPP

#include <vector>
#include <fstream>

#include "simulation.hpp"
#include "util.hpp"

struct Particle{
    long pid;
    double posX, posY, posZ;
    double hvX, hvY, hvZ;
    double velX, velY, velZ;
    double density;
    double accX, accY, accZ;

    Particle();

    Particle(std::ifstream& inputFile, int pid);

    int compute_grid_index(Simulation& sim);

    void write_particle_trace(std::ofstream& outputFile);

    void inc_part_dens (Particle & part_j, double const hSquared);

    void inc_part_acc (Particle &part_j, Simulation & sim, double const distanceSquared);

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
    Block (int bid, int blocks_x, int blocks_y, int blocks_z);
    int bid, index_i, index_j, index_k;
    std::vector<Particle> particles;
    std::vector<Block*> neighbours; // max size of 26

    // Methods 
    void reposition_particles();
    void compute_forces(); 
    void process_collisions();
    void move_particles(); 
    void process_boundaries();

    void block_part_col_xmin(Simulation &sim); // i = 0
    void block_part_col_xmax(Simulation &sim); // i = nx -1
    void block_part_col_ymin(Simulation &sim); // j = 0
    void block_part_col_ymax(Simulation &sim); // j = ny -1
    void block_part_col_zmin(Simulation &sim); // k = 0
    void block_part_col_zmax(Simulation &sim); // k = nz -1

    void boundint_xmin(Simulation &sim);
    void boundint_xmax(Simulation &sim);
    void boundint_ymin(Simulation &sim);
    void boundint_ymax(Simulation &sim);
    void boundint_zmin(Simulation &sim);
    void boundint_zmax(Simulation &sim);
};
#endif //COMP_ARCH_23_BLOCK_HPP
