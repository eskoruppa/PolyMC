#ifndef __WRITELOADRESTART_INCLUDED__
#define __WRITELOADRESTART_INCLUDED__

#include <string>    // string
#include <iostream>  // input output
#include <iomanip>   // std::setprecision
#include <stdexcept> // exceptions
#include <fstream>   // read file
#include <sstream>
#include <vector>
#include <armadillo> // armadillo
#include <algorithm>

#include "../Input/ReadFunctions.h"



struct RestartData;
struct RestartData {

    arma::cube      triads;
    arma::mat       pos;
    int             num_bp;
    std::string     type;
    std::string     sequence;
    long long int   snapshot;
    double          dLK;

};

std::vector<RestartData> loadrestart(std::string fn);




#endif
