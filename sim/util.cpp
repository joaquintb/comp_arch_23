#include "util.hpp"
#include "block.hpp"

// Function to show a .trz file
void display_trace(std::string filename) {
    using namespace std;

    struct traceParticle {
        int64_t id;
        double posX, posY, posZ;
        double hvX, hvY, hvZ;
        double velX, velY, velZ;
        double density;
        double accX, accY, accZ;
    };

    ifstream input_file(filename, ios::binary);

    if (!input_file.is_open()) {
        cerr << "Error opening the file." << '\n';
    }

    // Read the header (total number of blocks)
    int32_t totalBlocks;
    input_file.read(reinterpret_cast<char*>(&totalBlocks), sizeof(totalBlocks));

    cout << "Total number of blocks: " << totalBlocks << '\n';

    //cout << "Printing info about first 2 blocks: " << '\n';

    // Read information for each block
    for (int i = 0; i < totalBlocks; ++i) {
        // Read the number of particles in the block
        int64_t numParticles;
        input_file.read(reinterpret_cast<char*>(&numParticles), sizeof(numParticles));

        cout << "Block " << i + 1 << " - Number of particles: " << numParticles << '\n';

        // Read information for each particle in the block
        for (int j = 0; j < numParticles; ++j) {
            traceParticle particle;
            input_file.read(reinterpret_cast<char*>(&particle), sizeof(Particle));

            // Process the particle data as needed
            cout << "ID: " << particle.id << " | ";
            cout << "(" << particle.posX << ", " << particle.posY << ", " << particle.posZ << ")" << ", ";
            cout << "(" << particle.hvX << ", " << particle.hvY << ", " << particle.hvZ << ")" << ", ";
            cout << "(" << particle.velX << ", " << particle.velY << ", " << particle.velZ << ")" << ", ";
            cout << particle.density << ", ";
            cout << "(" << particle.accX << ", " << particle.accY << ", " << particle.accZ << ")" << '\n';
        }

    }

    input_file.close();
}