#ifndef __DUMP_XYZ_H__
#define __DUMP_XYZ_H__
#include "Dump.h"


class Dump_xyz;

class Dump_xyz: public Dump
{
protected:
    int counter;
    std::string center_conf;
    std::string rep;

public:
    Dump_xyz(Chain * ch, int N_dump, const std::string& filename,bool append=false, const std::string& center="COM", const std::string& representation="simple");
    ~Dump_xyz();

    void prod_dump();
    void final_dump();

protected:
    void xyz2file();

};




#endif
