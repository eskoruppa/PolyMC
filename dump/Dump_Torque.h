#ifndef __DUMP_TORQUE_H__
#define __DUMP_TORQUE_H__
#include "Dump.h"

////////////////////////////////////////////////////////////////////////////////
/*

*/

class Dump_Torque;

class Dump_Torque: public Dump
{
//protected:
//    int m_max;

public:
    Dump_Torque(Chain * ch, int N_dump, const std::string& filename, bool append=true);
    ~Dump_Torque();
    void prod_dump();
    void final_dump();
};

#endif
