#ifndef __EV_BEADS__
#define __EV_BEADS__

#include "../Chain.h"
#include "../SO3Methods.h"
#include "../ExtraFuncs.h"
#include "ProximityCluster.h"

#define _USE_MATH_DEFINES

#include <armadillo> // armadillo
#include <iostream>  // input output
#include <string>    // string
#include <algorithm>
#include <stdexcept> // exceptions
#include <fstream>   // read file
#include <cmath>     // math lib
#include <vector>    // vector
#include <cassert>   // assert

class ExcludedVolumeBead;

class ExcludedVolumeBead {

protected:

    int id;
    int id_refbead;
    double radius;
    int num_contained_bp;
    std::vector<int> contained_bp;
    arma::mat* all_pos;
    ProximityCluster* PC;

public:

    ExcludedVolumeBead(int EVBead_id, int refbead_id, int num_contained, double rad, arma::mat* pos);
    ~ExcludedVolumeBead();

    void set_ProximityCluster(ProximityCluster* bc);

    int get_id();
    int get_id_refbead();
    int get_num_contained_bp();
    std::vector<int>* get_contained_bp();

};

#endif
