#include "grid.hpp"

Grid::Grid(int size_x, int size_y, int size_z) {
    this->size_x = size_x;
    this->size_y = size_y;
    this->size_z = size_z;

    this->size = size_x * size_y * size_z;

    this->n_particles = 0;
    this->blocks = std::vector<Block>(this->size);
};

void Grid::populate(Simulation& sim, std::ifstream& inputFile) {
    int num_p = sim.get_num_p();

    for (int i = 0; i < num_p; i++) {
        Particle cur_particle = Particle(inputFile, i);
        int grid_index = cur_particle.compute_grid_index(sim);
        this->blocks[grid_index].particles.push_back(cur_particle);
    }

};

void Grid::display_grid(){
    for (auto block : this->blocks) {
        for (auto particle : block.particles) {
            std::cout << "ID: " << particle.pid  << " " << particle.posX << " " << particle.posY << " " << particle.posZ << std::endl;
        }
    }
}
