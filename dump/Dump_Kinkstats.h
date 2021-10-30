#ifndef __DUMP_KINKSTATS_H__
#define __DUMP_KINKSTATS_H__
#include "Dump.h"


class Dump_Kinkstats;

class Dump_Kinkstats: public Dump
{
protected:
    int option;
    int iID,fID;
    // 0 for exact
    // 1 for fuller
    // 2 for quick

//    double rkf_sum   = 0;
//    double rkf_count = 0;

    std::string fn_kinkpos;
    std::string fn_tancos;

    bool dump_fn_kinkpos = false;
    bool dump_fn_tancos  = false;

    std::vector<int> kinkpos;

    int num_omit_boundary;
    int first_considered;
    int last_considered;

    int neglect_below_kink_dist;

    int tancor_maxdist;
    arma::colvec tancor_sum;
    arma::colvec elastic_sum;
    arma::colvec tancor_count;


public:
    Dump_Kinkstats(Chain * ch, int N_dump, const std::string& filename, bool append=false);
    ~Dump_Kinkstats();

    void prod_dump();
    void final_dump();
};




#endif
