#ifndef __READFUNCTIONS_INCLUDED__
#define __READFUNCTIONS_INCLUDED__

#include <string>    // string
#include <iostream>  // input output
#include <stdexcept> // exceptions
#include <fstream>   // read file
#include <sstream>
#include <vector>
#include <armadillo> // armadillo
#include <algorithm>

void get_linelist(std::string &str, std::vector<std::string> &output);
void split_string(std::string &str, const std::string &delims, std::vector<std::string> &output);

#endif
