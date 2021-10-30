#ifndef __DUMP_TORQUETRAP_H__
#define __DUMP_TORQUETRAP_H__
#include "Dump.h"

////////////////////////////////////////////////////////////////////////////////
/*

*/

class Dump_TorqueTrap;

class Dump_TorqueTrap: public Dump
{

public:
    Dump_TorqueTrap(Chain * ch, int N_dump, const std::string& filename, bool append=true);
    ~Dump_TorqueTrap();
    void prod_dump();
    void final_dump();
};

#endif
