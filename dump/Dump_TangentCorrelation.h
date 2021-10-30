#ifndef __DUMP_TANCOR_H__
#define __DUMP_TANCOR_H__
#include "Dump.h"

////////////////////////////////////////////////////////////////////////////////
/*

*/

class Dump_TanCor;

class Dump_TanCor: public Dump
{
protected:
    arma::colvec tangent_dots;
    arma::colvec tangent_dots_count;
    int m_max;

public:
    Dump_TanCor(Chain * ch, int N_dump, const std::string& filename, int max_dist=50, bool append=true);
    ~Dump_TanCor();

    void prod_dump();
    void final_dump();
};

#endif
