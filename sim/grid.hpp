//
// Created by jtorresb on 5/11/23.
//

#ifndef COMP_ARCH_23_GRID_HPP
#define COMP_ARCH_23_GRID_HPP

#include <vector>
#include "block.hpp"

class Grid {
    public:
        Grid(int size_x, int size_y, int size_z);
    
    private:
        int size_x, size_y, size_z, n_particles;
        std::vector<Block> blocks;
};
#endif //COMP_ARCH_23_GRID_HPP
