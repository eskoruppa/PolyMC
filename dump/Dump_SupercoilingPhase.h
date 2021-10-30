#ifndef __DUMP_SUPERCOILINGPHASE_H__
#define __DUMP_SUPERCOILINGPHASE_H__
#include "Dump.h"


class Dump_SupercoilingPhase;

class Dump_SupercoilingPhase: public Dump
{
protected:
    double seg_size;

public:
    Dump_SupercoilingPhase(Chain * ch, int N_dump, const std::string& filename, bool append=false);
    ~Dump_SupercoilingPhase();

    void prod_dump();
    void final_dump();

protected:
    arma::colvec cal_writhe_dens(const arma::mat & writhemat);
    double       wr_locality(arma::colvec wrdens, double perc, double total);
    int          argmax(const arma::colvec & mat);
    arma::mat    neg2zero(const arma::mat& WM);
};




#endif
