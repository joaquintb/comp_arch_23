#ifndef COMP_ARCH_23_GRID_HPP
#define COMP_ARCH_23_GRID_HPP

#include "block.hpp"
#include "simulation.hpp"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <iostream>
#include <set>
#include <vector>

class Grid {
  public:
    Grid(int size_x, int size_y, int size_z);
    void populate(Simulation & sim, std::ifstream & inputFile);

    void set_neighbors();
    void increase_all_dens(Simulation & sim);
    void trans_all_dens(Simulation & sim);
    void increase_all_accs(Simulation & sim);
    void part_collisions();
    void motion();
    void part_box_collisions();
    void repos(Simulation & sim);
    void init_acc();
    void gen_output(std::ofstream & out);
    void simulate(int n_steps, Simulation & sim);

    // Trace comparison aux functions
    // bool cmp_trace(std::ifstream& trace);
    // void set_to_trace(std::ifstream& trace);

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

  private:
    int size_x, size_y, size_z, size;
    std::vector<Block> blocks;
};
#endif  // COMP_ARCH_23_GRID_HPP
