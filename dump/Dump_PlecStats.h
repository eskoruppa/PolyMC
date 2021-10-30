#ifndef __DUMP_PLECSTATS_H__
#define __DUMP_PLECSTATS_H__
#include "Dump.h"
#include "../PlecFinder/PlecFinder.h"

////////////////////////////////////////////////////////////////////////////////
/*
    Dump the extension of Writhe of a Tweezer Simulations
    Thus far only works in the constant linking number ensemble
*/

class Dump_PlecStats;
class Dump_PlecStats: public Dump
{
protected:

public:
    Dump_PlecStats(Chain * ch, int N_dump, const std::string& filename, bool append=true);
    ~Dump_PlecStats();

    void prod_dump();
    void final_dump();
};

#endif
