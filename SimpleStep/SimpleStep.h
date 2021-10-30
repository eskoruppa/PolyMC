#ifndef __SIMPLESTEP_H__
#define __SIMPLESTEP_H__

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
#include <random>    // random lib
#include <vector>    // STL vector




/*
    TODO:
        - allow moves to controle the move selection in a temperature dependent fashion.
          As of now temperature dependence is taken into account for the evaluation of the
          energies, but the proposed moves are not changed, which could lead to inefficient
          sampling.
*/

#define MCS_FORCE_ACTIVE            0
#define MCS_TORQUE_ACTIVE           1
#define MCS_TOPOLOGY_CLOSED         2
#define MCS_FIXED_LINK              3
#define MCS_FIXED_TERMINI           4
#define MCS_FIXED_TERMINI_RADIAL    5
#define MCS_FIXED_FIRST_ORIENTATION 6
#define MCS_FIXED_LAST_ORIENTATION  7

class MCStep;

class MCStep {

protected:
    string move_name="MCStep";

    // State Variables
    Chain * chain;
    vector<BPStep*> BPS;
    arma::mat*  pos;
    arma::cube* triads;


    int num_bp;
    int num_bps;

    double disc_len;

    // Statistics
    long long int count_step;
    long long int count_accept;

    // Random Devices
    std::random_device                      rd{};
    std::mt19937                            gen{rd()};
    std::normal_distribution<double>        normaldist{0.0,1.0};
    std::uniform_real_distribution<double>  uniformdist{0.0,1.0};
    std::uniform_int_distribution<>         randbinary{0, 1};
    /*
    Should those be static?
    */

    // moved intervals for excluded volume
    bool                requires_EV_check;
    vector<arma::ivec>  moved_intervals;
    bool                ev_active;
    long long int       count_ev_accept;
    arma::ivec          changed_bps;
    ExVol*              EV;

    arma::mat   avg_stiffmat;
    arma::mat   avg_covmat;

    bool suitablility_key[8] = {    false, // force_active
                                    false, // torque_active
                                    false, // topology_closed
                                    false, // link_fixed
                                    false, // fixed_termini
                                    false, // fixed_termini_radial
                                    false, // fixed_first_orientation
                                    false, // fixed_last_orientation
                                    };

public:
    MCStep(Chain * ch);
    ~MCStep();

// Monte Carlo Step (call main functionality)
public:
    // This method should not be overloaded
    bool MC();

/////////////////////////////////////
/*
    These functions need to be overloaded
*/
    virtual bool MC_move();
    void update_settings();
/////////////////////////////////////

// Get Moved Interval to take care of excluded volumes
public:
    void set_excluded_volume(ExVol* EVol);
    vector<arma::ivec>* get_moved_intervals();

// Check suitability of the move given the settings of Chain
    bool suitable();
    bool additional_criteria();


// Statistics
public:
    long long int get_steps();
    long long int get_accepts();
    long long int get_ev_accepts();
    double acceptance_rate();
    void print_acceptance_rate();



//protected:
//    arma::mat cal_avg_stiffmat();

};

#endif
