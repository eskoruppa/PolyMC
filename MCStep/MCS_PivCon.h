#ifndef __MCS_PIVCON_H__
#define __MCS_PIVCON_H__
#include "MCStep.h"

#include <cassert>
/*
Perhaps exclude the possibility to rotate te first triad.

*/

class MCS_PivCon;

class MCS_PivCon: public MCStep
{
protected:
    double    fac_sigma;
    arma::colvec  sigmas;
    int rge_id1,rge_id2;

    std::uniform_int_distribution<> genhinge1{0, 1};
    std::uniform_int_distribution<> genhingedist{0, 1};

public:
    MCS_PivCon(Chain * ch,const std::vector<long long int> & seedseq, int hingesize_min, int hingesize_max);
    MCS_PivCon(Chain * ch,const std::vector<long long int> & seedseq, int hingesize_min, int hingesize_max, int range_id1, int range_id2);
    ~MCS_PivCon();

    void update_settings();

// Monte Carlo Step (call main functionality)
public:
    bool MC_move();


};

#endif
