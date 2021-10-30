#ifndef __DUMP_TWIST_H__
#define __DUMP_TWIST_H__
#include "Dump.h"

#define DTW_STORE_TW 1000


class Dump_Twist;

class Dump_Twist: public Dump
{
protected:
    int iID,fID;
    arma::colvec store_tw;
    unsigned store_counter;

public:
    Dump_Twist(Chain * ch, int N_dump, const std::string& filename, bool append=true);
    Dump_Twist(Chain * ch, int N_dump, const std::string& filename, int start_id=0, int end_id=0, bool append=true);
    ~Dump_Twist();

    void init_store();
    void write2file();

    void prod_dump();
    void final_dump();
};

#endif
