#ifndef __CAL_LINKINGNUMBER_INCLUDED__
#define __CAL_LINKINGNUMBER_INCLUDED__

#include "../SO3Methods.h"
#include "../ExtraFuncs.h"
#include "../ExtraArma.h"

#include <armadillo> // armadillo
#include <iostream>  // input output
//#include <string>    // string
//#include <algorithm>
#include <stdexcept> // exceptions
//#include <fstream>   // read file
//#include <cmath>     // math lib
//#include <math.h>
//#include <vector>    // vector


#define LINKINKNUMBER_NORMZERO_THRESHOLD 1e-10


double fullerwrithe_vertex_angle    (const arma::colvec& a,const arma::colvec& b,const arma::colvec& c);
double fullerwrithe_triangle_area   (const arma::colvec& ez,const arma::colvec& t1,const arma::colvec& t2);
double fullerwrithe_writhesum       (const arma::mat& tans,const arma::colvec& z_dir);

double cal_closure_writhe(const arma::mat& pos,const arma::colvec& z_dir);

double cal_fullerwrithe             (const arma::mat& pos,const arma::colvec& z_dir);


//////////////////////////////////////////////////////////////////////////////
double       cal_langowskiwrithe(const arma::mat& pos,bool closed=false);

double       langowskiwrithe_quadrangle_area (const arma::colvec& r1,const arma::colvec& r2,const arma::colvec& r3,const arma::colvec& r4);
arma::colvec langowskiwrithe_nvec (const arma::colvec& a,const arma::colvec& b);
double       langowskiwrithe_angle_contribution (const arma::colvec& a,const arma::colvec& b);


double cal_langowskiwrithe_1b(const arma::mat& pos,bool closed=false);
double langowskiwrithe1b_Omega(const arma::colvec& r1,const arma::colvec& r2,const arma::colvec& r3,const arma::colvec& r4);
double langowskiwrithe1b_F(double t1,double t2,double a0,double cosbeta,double sinsqbeta);

#endif

