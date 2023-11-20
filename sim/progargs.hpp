#ifndef COMP_ARCH_23_PROGARGS_HPP
#define COMP_ARCH_23_PROGARGS_HPP

#include <algorithm>
#include <cctype>
#include <cmath>
#include <fstream>
#include <iostream>
#include <memory>
#include <numbers>
#include <span>
#include <string>
#include <typeinfo>
#include <vector>

int handle_num_args(int argc, std::string const & steps);
void check_files(std::ifstream & inputFile, std::string const & inputName,
                 std::ofstream & outputFile, std::string const & outputName);

#endif  // COMP_ARCH_23_PROGARGS_HPP
