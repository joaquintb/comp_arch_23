
Grid::Grid(int size_x, int size_y, int size_z) {
    this->size_x = size_x;
    this->size_y = size_y;
    this->size_z = size_z;
    this->n_particles = 0;
    this->blocks = std::vector<Block>();
}
