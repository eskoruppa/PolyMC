#ifndef __DUMP_ESENERGY__
#define __DUMP_ESENERGY__
#include "Dump.h"
#include "../Electrostatics/Electrostatics.h"

////////////////////////////////////////////////////////////////////////////////
/*
    Calculate Stiffness and Covariance Matrix
    both mean and individual base pair step stiffness (controlled by argument 'all_bps')
*/

class Dump_ESEnergy;

class Dump_ESEnergy: public Dump
{
protected:
//    bool electrostatics_active = false;
//    long long int counter;

    ElStat * elstat;
    bool kt_units = false;
    bool recal    = false;

public:
    Dump_ESEnergy(Chain * ch, int N_dump, const std::string& filename, ElStat * es, bool kt_units=false,bool recal=false, bool append=true);
    ~Dump_ESEnergy();
    void prod_dump();
};

#endif
