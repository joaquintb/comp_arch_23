//
// Created by jtorresb on 5/11/23.
//

#include "progargs.hpp"


void InputTest(int argc, char **argv){
    //handling alternative to direct use of char **argv
    std::span const args_view{argv, static_cast<std::size_t>(argc)};
    std::vector<std::string> const arguments{args_view.begin() + 1, args_view.end()};

    // ------------------------ READING INPUT FILE & ERROR HANDLING------------------------

    if(argc - 1 != 3){
        std::cerr << "Error: Invalid number of arguments: " << std::to_string(argc-1) << std::endl;
        std::exit(-1);
    }
    if(isdigit(std::to_integer((arguments[0]))) == false){
        std::cerr << "Error: time steps must be numeric." << std::endl;
        std::exit(-1);
    }
    if(std::to_integer(arguments[0]) < 0){
        std::cerr << "Error: Invalid number of time steps." << std::endl;
        std::exit(-2);
    }
    if (!arguments[1].is_open()) {
        std::cerr << "Error: Cannot open " <<arguments[1] <<  " for reading" << std::endl;
        std::exit(-3);
    }
    if (!arguments[2].is_open()) {
        std::cerr << "Error: Cannot open " <<arguments[2] <<  " for writing" << std::endl;
        std::exit(-4);
    }
}
