#ifndef __DUMP_MSD_H__
#define __DUMP_MSD_H__
#include "Dump.h"

////////////////////////////////////////////////////////////////////////////////
/*

*/

class Dump_MSD;

class Dump_MSD: public Dump
{
protected:
    arma::colvec sqdis;
    arma::colvec sqdis_count;
    int m_max;

public:
    Dump_MSD(Chain * ch, int N_dump, const std::string& filename, int max_dist=50, bool append=true);
    ~Dump_MSD();

    void prod_dump();
    void final_dump();
};

#endif
