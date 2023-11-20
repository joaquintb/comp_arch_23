#ifndef COMP_ARCH_23_BLOCK_HPP
#define COMP_ARCH_23_BLOCK_HPP

#include <vector>
#include <fstream>
#include <cmath>

#include "simulation.hpp"

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

struct Particle{
    long pid{};
    double posX{}, posY{}, posZ{};
    double hvX{}, hvY{}, hvZ{};
    double velX{}, velY{}, velZ{};
    double density{};
    double accX{}, accY{}, accZ{};

    Particle();

    Particle(std::ifstream& inputFile, int pid);

    int compute_grid_index(Simulation& sim) const;

    void write_particle_output(std::ofstream& outputFile) const;

    void inc_part_dens (Particle & part_j, double hSquared);

    void inc_part_acc (Particle &part_j, Simulation & sim, double distanceSquared);
};

class Block{
    // Attributes 
public:
    Block (int bid, int blocks_x, int blocks_y);
    int bid, index_i, index_j, index_k;
    std::vector<Particle> particles;
    std::vector<Block*> neighbours; // max size of 26

    void block_part_col_xmin(); // i = 0
    void block_part_col_xmax(); // i = nx -1
    void block_part_col_ymin(); // j = 0
    void block_part_col_ymax(); // j = ny -1
    void block_part_col_zmin(); // k = 0
    void block_part_col_zmax(); // k = nz -1

    void boundint_xmin();
    void boundint_xmax();
    void boundint_ymin();
    void boundint_ymax();
    void boundint_zmin();
    void boundint_zmax();
};
#endif //COMP_ARCH_23_BLOCK_HPP
