#ifndef __DUMP_STATE_H__
#define __DUMP_STATE_H__
#include "Dump.h"
#include "../ExcludedVolume/ExcludedVolume.h"


#include "../PlecFinder/PlecFinder.h"

/*
    TODO: Write to binary file rather than textfile.
*/


class Dump_State;

class Dump_State: public Dump
{
protected:
    ExVol * EV;
    bool EV_included;
    int counter;
    bool dump_pos = true;
    bool dump_triads;
    bool dump_Omegas;

public:
    Dump_State(Chain * ch, std::string mode, int N_dump, const std::string& filename,bool dump_triads=false,bool dump_Omegas=false, double EVrad=0, bool append=false);
    Dump_State(Chain * ch, std::string mode, int N_dump, const std::string& filename,bool dump_triads=false,bool dump_Omegas=false,bool append=false);
    Dump_State(Chain * ch, std::string mode, ExVol * EV, int N_dump, const std::string& filename,bool dump_triads=false,bool dump_Omegas=false,bool append=false);
    ~Dump_State();

    void prod_dump();
    void final_dump();
};




#endif
