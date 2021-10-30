#ifndef __INCLUDED_BRANCHBOXES__
#define __INCLUDED_BRANCHBOXES__

#include "../SO3Methods.h"
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


class BranchBox;

class BranchBox {
public:
    std::vector<int> xlim,ylim;
    arma::mat * WM;
    bool        active;
    double      writhe;
    bool        writhe_calculated;

public:
    BranchBox(const std::vector<int> & xlim,const std::vector<int> & ylim, arma::mat * WM);
    void set_WM(arma::mat * WM);
    void set_xlim(const std::vector<int> & xlim);
    void set_ylim(const std::vector<int> & ylim);
    void set_lims(const std::vector<int> & xlim,const std::vector<int> & ylim);
    double get_writhe();
    bool is_active();

protected:
    void check_finite_size();
};

#endif
