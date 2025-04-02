#ifndef __DUMP_CLOSUREHIST_H__
#define __DUMP_CLOSUREHIST_H__
#include "Dump.h"

////////////////////////////////////////////////////////////////////////////////
/*
    Calculate distance and angle histogram for Umbrella sampling if cyclization probabilities
*/

class Dump_ClosureHist;

class Dump_ClosureHist: public Dump
{
protected:
    bool active = true;
    double dist_lower;
    double dist_upper;
    double dist_range;

    bool bin_costheta;
    double angle_upper;
    double angle_lower;
    double angle_range;

    double twist_upper;
    double twist_lower;
    double twist_range;

    bool dist_active  = false;
    bool angle_active = false;
    bool twist_active = false;
    int num_active = 0;
    int combination_id = -1;
    int num_dist_bins;
    int num_angle_bins;
    int num_twist_bins;

    int num_bins_1 = 0;
    int num_bins_2 = 0;
    int num_bins_3 = 0;
    
    arma::ivec hist;
    arma::imat mathist;
    arma::icube cubehist;
    long long int invalid_counter;


public:
    Dump_ClosureHist(
        Chain * ch, 
        int N_dump, 
        const std::string& filename, 
        int num_dist_bins, 
        int num_angle_bins, 
        int num_twist_bins, 
        double dist_lower, 
        double dist_upper = -1, 
        double angle_lower = 0, 
        double angle_upper = 0,
        double twist_lower = 0,
        double twist_upper = 0,
        bool bin_costheta = false);
    ~Dump_ClosureHist();

    void prod_dump();
    void final_dump();

    bool is_active();
};

#endif
