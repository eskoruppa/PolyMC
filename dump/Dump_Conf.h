#ifndef __DUMP_CONF_H__
#define __DUMP_CONF_H__
#include "Dump.h"


class Dump_Conf;

class Dump_Conf: public Dump
{

public:
    Dump_Conf(Chain * ch, int N_dump, const std::string& filename,bool append=false);
    ~Dump_Conf();

    void prod_dump();
    void final_dump();
};




#endif
