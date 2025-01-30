#ifndef __DUMP_CLOSURE_H__
#define __DUMP_CLOSURE_H__
#include "Dump.h"

////////////////////////////////////////////////////////////////////////////////
/*
    Calculate distance and angle histogram and counts occurances of closure
*/

#define CLOSURE_NUM_STORES 50

class Dump_Closure;

class Dump_Closure: public Dump
{
protected:

    bool use_dist_criterion  = true;
    bool use_angle_criterion = true;
    bool use_twist_criterion = true;

    double dist_threshold;
    double angle_threshold;
    double twist_threshold;

    bool active;

    int store_counter;
    arma::mat stored_dumps;

public:
    Dump_Closure(
        Chain * ch, 
        int N_dump, 
        const std::string& filename, 
        double dist_threshold,
        double angle_threshold,
        double twist_threshold
        );
    ~Dump_Closure();

    void prod_dump();
    void final_dump();

    void init_store();
    void write2file();
    bool is_active();
};

#endif
