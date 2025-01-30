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

    // bool dist_active = true;
    bool angle_active = false;
    int num_dist_bins;
    int num_angle_bins;
    
    arma::ivec hist;
    // arma::colvec angles;
    arma::imat mathist;
    long long int invalid_counter;


public:
    Dump_ClosureHist(
        Chain * ch, 
        int N_dump, 
        const std::string& filename, 
        int num_dist_bins, 
        int num_angle_bins, 
        double dist_lower, 
        double dist_upper = -1, 
        double angle_lower = 0, 
        double angle_upper = 0,
        bool bin_costheta = true);
    ~Dump_ClosureHist();

    void prod_dump();
    void final_dump();

    bool is_active();
};

#endif
