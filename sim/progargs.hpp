#ifndef COMP_ARCH_23_PROGARGS_HPP
#define COMP_ARCH_23_PROGARGS_HPP

#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
#include <memory>
#include <algorithm>
#include <numbers>
#include <span>
#include <typeinfo>
#include <cctype>
#include <string>

int handle_num_args(int argc, std::string const & steps);
void check_files(std::ifstream & inputFile, std::string const & inputName,
                 std::ofstream & outputFile, std::string const & outputName);

#endif //COMP_ARCH_23_PROGARGS_HPP
