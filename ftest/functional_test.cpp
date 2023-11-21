#include <gtest/gtest.h>

#include "../sim/block.hpp"
#include "../sim/grid.hpp"
#include "../sim/simulation.hpp"
#include "../sim/progargs.hpp"
#include <cmath>
#include <fstream>
#include <iostream>
#include <string>

std::string ReadFileToString(const std::string& file_path) {
    const std::ifstream file(file_path);
    std::stringstream buffer;
    buffer << file.rdbuf();
    return buffer.str();
}

TEST(SimulationFunctionalTest, SmallInputFunctionalTest) {
    const int num_iter = 5;
    for (int i = 1; i <= num_iter; i++) {
        const std::string input_file = "small.fld";
        const std::string output_file = "out_small_ftest.fld";

        std::ifstream input_stream(input_file, std::ios::binary);
        std::ofstream output_stream(output_file, std::ios::binary | std::ios::trunc);
        const int n_steps = i;

        auto ppm = read_binary_value<float>(input_stream);
        auto num_p = read_binary_value<int>(input_stream);

        Simulation sim = Simulation(ppm, num_p);
        Grid grid = Grid(sim.n_x, sim.n_y, sim.n_z);

        grid.populate(sim, input_stream);
        grid.set_neighbors();
        input_stream.close();
        grid.simulate(n_steps, sim);

        write_binary_value(static_cast<float>(ppm), output_stream);
        write_binary_value(static_cast<int>(num_p), output_stream);
        grid.gen_output(output_stream);
        output_stream.close();

        const std::string reference_file = "out/small-" + std::to_string(i) + ".fld";

        const std::string output_content = ReadFileToString(output_file);
        const std::string reference_content = ReadFileToString(reference_file);
        EXPECT_EQ(output_content, reference_content) << "The contents of " << output_file 
                                                     << " and " << reference_file 
                                                     << " do not match.";
        
        std::cout << "Passed test " << i << " of 5 for SmallInputFunctionalTest." << "\n"; 
    }
}

TEST(SimulationFunctionalTest, LargeInputFunctionalTest) {
    const int num_iter = 5;
    for (int i = 1; i <= num_iter; i++) {
        const std::string input_file = "large.fld";
        const std::string output_file = "out_large_ftest.fld";

        std::ifstream input_stream(input_file, std::ios::binary);
        std::ofstream output_stream(output_file, std::ios::binary | std::ios::trunc);
        const int n_steps = i;

        auto ppm = read_binary_value<float>(input_stream);
        auto num_p = read_binary_value<int>(input_stream);

        Simulation sim = Simulation(ppm, num_p);
        Grid grid = Grid(sim.n_x, sim.n_y, sim.n_z);

        grid.populate(sim, input_stream);
        grid.set_neighbors();
        input_stream.close();
        grid.simulate(n_steps, sim);

        write_binary_value(static_cast<float>(ppm), output_stream);
        write_binary_value(static_cast<int>(num_p), output_stream);
        grid.gen_output(output_stream);
        output_stream.close();

        const std::string reference_file = "out/large-" + std::to_string(i) + ".fld";

        const std::string output_content = ReadFileToString(output_file);
        const std::string reference_content = ReadFileToString(reference_file);
        EXPECT_EQ(output_content, reference_content) << "The contents of " << output_file 
                                                     << " and " << reference_file 
                                                     << " do not match.";
        
        std::cout << "Passed test " << i << " of 5 for LargeInputFunctionalTest." << "\n"; 
    }
}
