#ifndef __DUMP_TS_H__
#define __DUMP_TS_H__
#include "Dump.h"

////////////////////////////////////////////////////////////////////////////////
/*

*/
#define ETS_UNDEFINED -1
#define ETS_FULLER    1
#define ETS_LANGOWSKI 2
#define ETS_BOTH      3

class Dump_EffTorsStiff;

class Dump_EffTorsStiff: public Dump
{
protected:
    double dLK, dLK_sq, Wr, Wr_sq, Tw, Tw_sq;
    double Wr_ful,Wr_sq_ful;
    long long int counter;
    int     iID,fID;
    std::string  method;
    int     method_id;
    arma::colvec dir;

public:
    Dump_EffTorsStiff(Chain * ch, int N_dump, const std::string& filename, const arma::colvec& force_dir, const std::string& writhe_method="fuller", int start_id=0, int end_id=0, bool append=true);
    ~Dump_EffTorsStiff();

    void prod_dump();
    void final_dump();
};

#endif
