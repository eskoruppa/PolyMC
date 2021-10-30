#ifndef __DUMP_DISTMAP_H__
#define __DUMP_DISTMAP_H__
#include "Dump.h"


class Dump_DistMap;

class Dump_DistMap: public Dump
{
protected:
    double density;

public:
    Dump_DistMap(Chain * ch, int N_dump, const std::string& filename, double seg_size, bool append=false);
    ~Dump_DistMap();

    void prod_dump();
    void final_dump();
};




#endif
