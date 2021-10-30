#ifndef __INCLUDED_PROXMITYPLECFINDER__
#define __INCLUDED_PROXMITYPLECFINDER__

#include "../Chain.h"
#include "../SO3Methods.h"
#include "../ExtraFuncs.h"
#include "../ExtraArma.h"
#include "BranchBoxes.h"

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
#include <tuple>     // for multiple function returns

//#define PLECFINDER_DEBUG


class ProximityPlecFinder;

class ProximityPlecFinder {
protected:
    int     AUTODETECT_ITERATIONS = 100;

    double threshold_dist        = 0;
    int    min_seg_dist          = 0;
    double frac_closest          = 0.9;

    bool   autodetect_thresholds = false;
    int    autodetect_count      = 0;
    std::vector<double> autodetected_thresholds;

    arma::ivec closest_id;
    arma::vec  closest_dist;


public:

    ProximityPlecFinder();
    ProximityPlecFinder(double threshold_dist,int min_seg_dist,double frac_closest);
    ~ProximityPlecFinder();

    double get_xp(const arma::mat * pos);
    std::vector<int> find_plecs(const arma::mat * pos);
    std::vector<int> find_plecs_Neuman(const arma::mat * pos);

private:

    void cal_closest(const arma::mat * pos);

};

#endif
