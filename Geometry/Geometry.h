#ifndef __GEOMETRY__
#define __GEOMETRY__

#include <armadillo>
#include <cmath>     // math lib

double distance_point_lineseg(const arma::colvec p, const arma::colvec q1, const arma::colvec q2);
double distance_point_semi_infinite_lineseg(const arma::colvec p, const arma::colvec q1, const arma::colvec v);

double distance_parallel_linesegs_finite_semiinf(const arma::colvec p1, const arma::colvec p2, const arma::colvec q, const arma::colvec v);

bool vectors_parallel(const arma::colvec v1, const arma::colvec v2, double tol=1e-10);

#endif
