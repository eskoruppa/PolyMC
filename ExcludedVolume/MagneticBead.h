#ifndef __EXCLUDED_VOLUME__
#define __EXCLUDED_VOLUME__

#include "../Chain.h"
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


/*
    DEBUG FLAG
*/
//#define DEBUG_MAGNETICBEAD


class MagneticBead;

class MagneticBead {

protected:

    Chain * chain;
    std::vector<BPStep*> BPS;
    arma::mat*  bp_pos;
    arma::cube* triads;

    arma::mat*  bp_pos_backup;
    arma::cube* triads_backup;

    int num_bp;
    int num_bps;
    double disc_len;
    bool closed_topology;

    double      EV_rad;
    double      bead_rad;
    double      min_dist;

    int         attach_bp;
    arma::vec   rel_vec;

    int counter;
    int counter_reject;




public:
    MagneticBead(Chain * ch, arma::mat * bp_pos_backup, arma::cube * triads_backup, double bead_radius, double EV_radius, int attachment_bp ,const arma::vec & vec_triad2bead);
    ~MagneticBead();

public:
    bool check_bead_EV_double();
    bool check_bead_EV_current();

    arma::vec get_bead_pos();
    arma::vec get_bead_pos_backup();

public:
    double rejection_rate();


protected:
    double doubleMove(const arma::vec & beadpos_past, const arma::vec & beadpos_curr, int bp_id);

};
#endif
