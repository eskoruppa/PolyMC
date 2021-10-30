#ifndef __DUMP_WRITHEMAP_H__
#define __DUMP_WRITHEMAP_H__
#include "Dump.h"


class Dump_WritheMap;

class Dump_WritheMap: public Dump
{
protected:
    double seg_size;

public:
    Dump_WritheMap(Chain * ch, int N_dump, const std::string& filename, double seg_size, bool append=false);
    ~Dump_WritheMap();

    void prod_dump();
    void final_dump();
};




#endif
