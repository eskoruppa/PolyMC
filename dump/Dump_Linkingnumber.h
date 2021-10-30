#ifndef __DUMP_LINKINGNUMBER_H__
#define __DUMP_LINKINGNUMBER_H__
#include "Dump.h"

#include "../Calcs/LinkingNumber.h"


////////////////////////////////////////////////////////////////////////////////
/*
    options determined how the writhe is calculated:
    - "exact"       calculates the exact writhe with the langowski method (1a)
    - "fuller"      uses fuller writhe that is valid modulo 1
    - "quick"       If the linking number is a conserved quantity in this simulation, the writhe
                    will be calculated as the difference between the excess linking number and excess twist.
*/

class Dump_Linkingnumber;

class Dump_Linkingnumber: public Dump
{
protected:
    int option;
    int iID,fID;
    // 0 for exact
    // 1 for fuller
    // 2 for quick

public:
    Dump_Linkingnumber(Chain * ch, int N_dump, const std::string& filename, const std::string& options="exact", bool append=false);
    Dump_Linkingnumber(Chain * ch, int N_dump, const std::string& filename, int start_id, int end_id, const std::string& options="exact", bool append=false);
    ~Dump_Linkingnumber();

    void prod_dump();
};




#endif
