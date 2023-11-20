#include <iostream>
#include <fstream>

struct Header {
    float ppm;
    int num_particles;
};

void decreaseNumParticles(const char* moreFileName) {
    // Step 1: Open more.fld in binary mode
    std::fstream moreFile(moreFileName, std::ios::binary | std::ios::in | std::ios::out);

    // Step 2: Decrease the num_particles value by 1
    Header currentHeader;
    moreFile.read(reinterpret_cast<char*>(&currentHeader), sizeof(Header));

    // Assuming you want to decrease the num_particles by 1
    currentHeader.num_particles++;

    // Step 3: Seek to the beginning of the file
    moreFile.seekp(0);

    // Step 4: Write the updated header back to the file
    moreFile.write(reinterpret_cast<char*>(&currentHeader), sizeof(Header));

    // File will be closed automatically when moreFile goes out of scope
}

int main() {
    const char* moreFileName = "less.fld";

    decreaseNumParticles(moreFileName);

    std::cout << "num_particles increased successfully." << std::endl;

    return 0; //
}
