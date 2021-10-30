#ifndef __DUMP_E2ED_H__
#define __DUMP_E2ED_H__
#include "Dump.h"

////////////////////////////////////////////////////////////////////////////////

class Dump_EndToEndDistance;

class Dump_EndToEndDistance: public Dump
{
public:
    Dump_EndToEndDistance(Chain * ch, int N_dump, const std::string& filename, bool append=true);
    ~Dump_EndToEndDistance();

    void prod_dump();
    void final_dump();
};

#endif
