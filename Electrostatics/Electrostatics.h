#ifndef __EXCLUDED_ELECTROSTATICS__
#define __EXCLUDED_ELECTROSTATICS__

#include "../Chain.h"
#include "ESPotential.h"
#include "../SO3Methods.h"
#include "../ExtraFuncs.h"
#include "../ExcludedVolume/ExcludedVolume.h"
#include "../Input/InputRead.h"

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
    Algorithm Flags
*/
#define USE_STEPSKIP
#define RECAL_FULL_EVERY 1000
#define MAX_ENERGY_DIFF  1e-10

/*
    EV Type defines are given in ExcludedVolume.h
*/




/*
    DEBUG FLAG
*/
//#define DEBUG_ELSTAT


class ElStat;

class ElStat {

protected:

    Chain * chain;
    std::vector<BPStep*> BPS;
    arma::mat*  bp_pos;
    arma::cube* triads;

    ESPotential * espot;

    int     num_bp;
    int     num_bps;
    double  disc_len;
    double  kT;


    bool closed_topology;

    double rho_max;
    double rho_min;

    int num_midpoint;
    int half_num_midpoint;

    arma::mat  midpoint_pos;
    arma::mat  midpoint_pos_backup;

    arma::mat  pair_interactions;
    arma::mat  pair_interactions_backup;

    double  neighbor_skip_dist;
    int     neighbor_skip_num;

    bool pending_check;
    long long counter;


public:

    ElStat(Chain * ch, ESPotential * espot, double rho_max, double rho_min, double neighbor_skip_dist);
    ~ElStat();

public:

    void    eval_full();
    double  cal_beta_dE(const std::vector<arma::ivec>* moved);
//    void    set_backup(arma::mat * midpoint_pos, arma::mat * pair_interactions);

    double get_current_energy(bool kt_units=false,bool recal=false);

    arma::mat* get_midpoint_pos();
    arma::mat* get_midpoint_pos_backup();
    arma::mat* get_pair_interactions();
    arma::mat* get_pair_interactions_backup();



protected:
//    void    cal_midpoint_pos(const std::vector<arma::ivec>* moved);
    void    cal_midpoint_pos(const std::vector<arma::ivec> & midpoint_intervals);
    void    cal_midpoint_pos();

    void    cal_midpoint_intervals(    const std::vector<arma::ivec>* moved,
                                       std::vector<arma::ivec>& midpoint_intervals);

    int     impose_skip(int a, int b);
    bool    in_neighborskip_distance(int a, int b);

protected:
    double eval_dE_inter_interval(int a1, int a2, int b1, int b2);
    double eval_dE_intra_interval(int a1, int a2);
    double eval_dE_bead_on_interval(int a, int b1, int b2);

    double eval_pair(int a, int b);
    double eval_pair(int a, int b, double dist);

public:
// Backup Handling
    void    revert_to_backup();
    void    set_current_as_backup(bool recal=false);

protected:
// check energy
    void recal_full_energy();

// debug functions
public:
    double  check_eval_pair(int a, int b);
    bool    check_midpoint_pos();
    bool    check_pair_interactions();

};


#endif
