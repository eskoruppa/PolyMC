#ifndef __DUMP_STIFF_H__
#define __DUMP_STIFF_H__
#include "Dump.h"

////////////////////////////////////////////////////////////////////////////////
/*
    Calculate Stiffness and Covariance Matrix
    both mean and individual base pair step stiffness (controlled by argument 'all_bps')
*/

class Dump_Stiff;

class Dump_Stiff: public Dump
{
protected:
    int  counter;
    bool all;
    arma::cube covs;
    arma::mat  cov;

public:
    Dump_Stiff(Chain * ch, int N_dump, const std::string& filename, bool all_bps=false);
    ~Dump_Stiff();

    void prod_dump();
    void final_dump();
};

#endif
