#ifndef __DUMP_HELICITY_H__
#define __DUMP_HELICITY_H__
#include "Dump.h"


////////////////////////////////////////////////////////////////////////////////
/*
    options determined how the writhe is calculated:
    - "exact"       calculates the exact writhe with the langowski method (1a)
    - "fuller"      uses fuller writhe that is valid modulo 1
    - "quick"       If the linking number is a conserved quantity in this simulation, the writhe
                    will be calculated as the difference between the excess linking number and excess twist.
*/

class Dump_Helicity;

class Dump_Helicity: public Dump
{
protected:
    int    N_dump2file;
    int    dump2file_counter;
    double hel_sum;


public:
    Dump_Helicity(Chain * ch, int N_dump, int N_dump2file, const std::string& filename, bool append=false);
    ~Dump_Helicity();

    void prod_dump();
};




#endif
