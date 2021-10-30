#ifndef __DUMP_ENERGY__
#define __DUMP_ENERGY__
#include "Dump.h"

////////////////////////////////////////////////////////////////////////////////
/*
    Calculate Stiffness and Covariance Matrix
    both mean and individual base pair step stiffness (controlled by argument 'all_bps')
*/

class Dump_Energy;

class Dump_Energy: public Dump
{
protected:
    double z,z_sq;
    long long int counter;
    int iID,fID;
    bool stepwise;

public:
    Dump_Energy(Chain * ch, int N_dump, const std::string& filename, bool stepwise=true, bool append=true);
    ~Dump_Energy();

    void prod_dump();
};

#endif
