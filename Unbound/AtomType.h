#ifndef __ATOMTYPE_INCLUDED__
#define __ATOMTYPE_INCLUDED__

#include "../ExtraFuncs.h"

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


struct AtomType;
struct AtomType {
    int    id;
    double hard_repulsion_cuttoff;
    std::vector<double> params;
    bool locked = false;
};

struct TypeGroup;
struct TypeGroup {
    AtomType            type;
    int                 id;
    std::vector<int>    elements;

    bool                involved_in_pairinteraction = false;

//    bool                trial_pending = false;
    // New elements suggested by element switch move
    std::vector<int>    elements_trial;
    std::vector<int>    elements_removed;
    std::vector<int>    elements_added;
    std::vector<int>    elements_remaining;

    // Invervals for Regular MC moves
    std::vector<int>    unmoved;
    std::vector<int>    moved;
    std::vector<int>    individual;

    int                 num_moved=0;

//    void trial_add(const std::vector<int> & add);
//    void trial_remove(const std::vector<int> & remove);

};


#endif
