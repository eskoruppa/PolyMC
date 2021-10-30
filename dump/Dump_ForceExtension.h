#ifndef __DUMP_FE_H__
#define __DUMP_FE_H__
#include "Dump.h"

////////////////////////////////////////////////////////////////////////////////


class Dump_ForceExtension;

class Dump_ForceExtension: public Dump
{
protected:
    double z,z_sq;
    long long int counter;
    int iID,fID;

public:
    Dump_ForceExtension(Chain * ch, int N_dump, const std::string& filename, int start_id=0, int end_id=0, bool append=true);
    ~Dump_ForceExtension();

    void prod_dump();
    void final_dump();
};

#endif
