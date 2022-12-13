#ifndef __CLOSURECROSSINGS__
#define __CLOSURECROSSINGS__

#define _USE_MATH_DEFINES
#include <armadillo> // armadillo

#include "../Geometry/Geometry.h"

double front_closure(   const arma::colvec p1,
                        const arma::colvec p2,
                        const arma::colvec r0,
                        const arma::colvec y,
                        double tol=1e-12);
double back_closure(    const arma::colvec p1,
                        const arma::colvec p2,
                        const arma::colvec r_l1,
                        const arma::colvec r_l2,
                        const arma::colvec y,
                        double tol=1e-12);


#endif
