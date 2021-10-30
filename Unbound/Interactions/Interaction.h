#ifndef __INTERACTION_INCLUDED__
#define __INTERACTION_INCLUDED__

#include "../../Chain.h"
#include "../../SO3Methods.h"
#include "../../ExtraFuncs.h"
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


class Interaction;
class Interaction {

protected:
    Chain *      chain;
    arma::mat**  new_pos;
    arma::cube** new_trds;
    arma::mat**  old_pos;
    arma::cube** old_trds;


public:

    Interaction(Chain * ch, arma::mat** new_pos, arma::cube** new_trds, arma::mat**  old_pos, arma::cube** old_trds);
    ~Interaction();

    double EvalEnergy(int id1, int id2, bool check_old, bool check_new);

    virtual double CalEnergy();





};

#endif
