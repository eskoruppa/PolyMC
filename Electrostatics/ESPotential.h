#ifndef __EXCLUDED_ESPOTENTIAL__
#define __EXCLUDED_ESPOTENTIAL__

#include "../Chain.h"
#include "../SO3Methods.h"
#include "../ExtraFuncs.h"
#include "../ExcludedVolume/ExcludedVolume.h"

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

/*
    CONSTANTS
*/
#define CHARGE_PER_BP           -2
#define DISTANCE_BETWEEN_BP     0.34

#define HALF_PI                 1.5707963267948966192313216916397514420985846996875529104874722961
#define ZERO_CROSS              1e-10

/*
    DEBUG FLAG
*/
//#define DEBUG_ESPOTENTIAL


class ESPotential;

class ESPotential {

protected:

    std::string type = "ESPotential";

    Chain * chain;
    double  disc_len;
    double  kT;

    std::vector<double> & params;

    double  integral_dx;
    double  rho_max;
    int     num_discretizations;

    double  charge_density;
    double  prefactor;

    bool    tabulate;
    bool    use_costheta;
    int     table_N;
    int     table_intervals;
    double  min_dist;
    double  max_dist;
    double  memory_limit;

    std::vector<arma::cube> table;
    double          table_dist_step;
    double          table_dist_range;
    arma::colvec    table_dist_vals;

    double          table_costheta_step;
    arma::colvec    table_costheta_vals;


public:

    ESPotential(Chain * ch, std::vector<double> & params, double integral_dx, double rho_max);
    ~ESPotential();

public:
    bool   tabulate_interactions(int table_N,double min_dist,double max_dist,double memory_limit,bool use_costheta=true);

    std::string get_type();
    bool        tabulating();

public:

    double  segmentpair_energy( const arma::colvec & r,
                                const arma::colvec & t1,
                                const arma::colvec & t2);

    double  segmentpair_energy( const arma::colvec & r,
                                const arma::colvec & t1,
                                const arma::colvec & t2,
                                const double dist);

protected:
    double  eval_double_integral(   const arma::colvec & r,
                                    const arma::colvec & t1,
                                    const arma::colvec & t2,
                                    const double dist);

    double  eval_double_integral(   const arma::colvec & r,
                                    const arma::colvec & t1,
                                    const arma::colvec & t2,
                                    const double dist,
                                    double dx);
    virtual double set_prefactor();
    virtual double integrant(double r);


protected:
    void init_tabulation_costheta();
    void init_tabulation_theta();

    double  eval_table( const arma::colvec & r,
                        const arma::colvec & t1,
                        const arma::colvec & t2,
                        const double dist);
    double  eval_table_theta(   const arma::colvec & r,
                                const arma::colvec & t1,
                                const arma::colvec & t2,
                                const double dist);


    inline double   id2theta(int id);
    inline int      costheta2id(double costheta);

    inline void     set_table_val(double val, int d, int t1, int t2, int t3);
    inline double   get_table_val(            int d, int t1, int t2, int t3);

};

inline double   ESPotential::id2theta(int id) {
    return std::acos( -1+id*table_costheta_step );
}

inline int      ESPotential::costheta2id(double costheta) {

}

inline void ESPotential::set_table_val(double val, int d, int t1, int t2, int t3) {
    table[d](t1,t2,t3) = val;
}

inline double ESPotential::get_table_val(int d, int t1, int t2, int t3) {
    return table[d](t1,t2,t3);
}






#endif
