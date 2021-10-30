#ifndef __DUMP_LB_H__
#define __DUMP_LB_H__
#include "Dump.h"

////////////////////////////////////////////////////////////////////////////////
/*
    Calculate Stiffness and Covariance Matrix
    both mean and individual base pair step stiffness (controlled by argument 'all_bps')
*/

class Dump_PersLen;

class Dump_PersLen: public Dump
{
protected:
    arma::colvec tangent_dots;
    arma::colvec tangent_dots_count;
    int m_max;
    bool cal_torsional_perslen = false;
    arma::colvec cosom3_sums;
    arma::colvec cosom3_sums_count;

public:
    Dump_PersLen(Chain * ch, int N_dump, const std::string& filename, int max_dist=50,bool cal_tors=false,bool append=true);
    ~Dump_PersLen();

    void prod_dump();
    void final_dump();
};

#endif
