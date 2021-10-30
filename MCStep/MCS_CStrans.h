#ifndef __MCS_CSTRANS_H__
#define __MCS_CSTRANS_H__
#include "MCStep.h"

#include <cassert>
/*
Perhaps exclude the possibility to rotate te first triad.

*/

class MCS_CStrans;

class MCS_CStrans: public MCStep
{
protected:
    arma::mat covmat;
    double    fac;

    double       sigma;
    arma::colvec trans_dir;
    bool         fixed_dir;

    std::uniform_int_distribution<> genhinge1{0, 1};
    std::uniform_int_distribution<> genhingedist{0, 1};

public:
    MCS_CStrans( Chain * ch,const std::vector<long long int> & seedseq, int selrange_min, int selrange_max, arma::colvec dir={0,0,0} );
    ~MCS_CStrans();

    void update_settings();

// Monte Carlo Step (call main functionality)
public:
    bool MC_move();


};

#endif
