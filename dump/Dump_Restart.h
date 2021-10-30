#ifndef __DUMP_RESTART_H__
#define __DUMP_RESTART_H__
#include "Dump.h"

/*
    TODO: Write to binary file rather than textfile.
*/


class Dump_Restart;

class Dump_Restart: public Dump
{
protected:
    int counter;
    bool last;

public:
    Dump_Restart(Chain * ch, int N_dump, const std::string& filename, bool append=false, bool last_snapshot=false);
    ~Dump_Restart();

    void prod_dump();
    void final_dump();
};


#endif



