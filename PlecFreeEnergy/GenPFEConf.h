#ifndef __INCLUDED_GENPFECONF__
#define __INCLUDED_GENPFECONF__

#include "../SO3Methods.h"
#include "../ExtraFuncs.h"
#include "../ExtraArma.h"

#define _USE_MATH_DEFINES

#include <armadillo> // armadillo
#include <iostream>  // input output
#include <algorithm>
#include <stdexcept> // exceptions
#include <fstream>   // read file
#include <cmath>     // math lib
#include <vector>    // vector


int gen_PFE_conf(arma::mat * pos, arma::cube * triads, int num_segs, double  disc_len, double termini_dist, double EV_rad=0, double loop_frac=2./3);



#endif
