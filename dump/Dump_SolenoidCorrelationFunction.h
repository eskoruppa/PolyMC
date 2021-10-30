#ifndef __DUMP_SOLCOR_H__
#define __DUMP_SOLCOR_H__
#include "Dump.h"

////////////////////////////////////////////////////////////////////////////////
/*

*/

class Dump_SolCor;

class Dump_SolCor: public Dump
{
protected:
    arma::mat cors;
    arma::mat cors_count;
    int nmax,kmax;
    arma::colvec zdir;
    bool use_tangents=true;

public:
    Dump_SolCor(Chain * ch, int N_dump, const std::string& filename, arma::colvec z_dir,int k_max, int n_max=50, bool append=true);
    ~Dump_SolCor();

    void prod_dump();
    void final_dump();
};

#endif
