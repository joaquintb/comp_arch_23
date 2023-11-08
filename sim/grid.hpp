#ifndef COMP_ARCH_23_GRID_HPP
#define COMP_ARCH_23_GRID_HPP

#include <vector>
#include "block.hpp"
#include "simulation.hpp"

class Grid {
    public:
        Grid(int size_x, int size_y, int size_z);
        void populate(Simulation& sim, std::ifstream& inputFile);
    
    private:
        int size_x, size_y, size_z, n_particles, size;
        std::vector<Block> blocks;
};
#endif //COMP_ARCH_23_GRID_HPP
