#include "progargs.hpp"


void inputTest(int argc, char **argv){
    //handling alternative to direct use of char **argv
    std::span const args_view{argv, static_cast<std::size_t>(argc)};
    std::vector<std::string> const arguments{args_view.begin() + 1, args_view.end()};

    std::ifstream input (arguments[1]);
    std::ifstream output (arguments[2]);

    bool numeric = true;

    // ------------------------ READING INPUT FILE & ERROR HANDLING------------------------

    if(argc - 1 != 3){
        std::cerr << "Error: Invalid number of arguments: " << std::to_string(argc-1) << std::endl;
        std::exit(-1);
    }
    if(numeric == false){
        std::cerr << "Error: time steps must be numeric." << std::endl;
        std::exit(-1);
    }
    if(std::stoi(arguments[0]) < 0){
        std::cerr << "Error: Invalid number of time steps." << std::endl;
        std::exit(-2);
    }
    if (!input.is_open()) {
        std::cerr << "Error: Cannot open " <<arguments[1] <<  " for reading" << std::endl;
        std::exit(-3);
    }
    if (!output.is_open()) {
        std::cerr << "Error: Cannot open " <<arguments[2] <<  " for writing" << std::endl;
        std::exit(-4);
    }
}