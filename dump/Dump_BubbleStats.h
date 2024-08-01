#ifndef __DUMP_BUBBLESTATS_H__
#define __DUMP_BUBBLESTATS_H__
#include "Dump.h"

#define DBUB_STORE_NUM 1000

////////////////////////////////////////////////////////////////////////////////

class Dump_BubbleStats;

class Dump_BubbleStats: public Dump
{
protected:
    // int iID,fID;
    arma::colvec store_num_bubsegs;
    arma::colvec store_num_bubs;
    unsigned store_counter;

public:
    Dump_BubbleStats(Chain * ch, int N_dump, const std::string& filename, bool append=true);
    // Dump_BubbleStats(Chain * ch, int N_dump, const std::string& filename, int start_id, int end_id, bool append=true);
    ~Dump_BubbleStats();

    void init_store();
    void write2file();

    void prod_dump();
    void final_dump();
};

#endif
