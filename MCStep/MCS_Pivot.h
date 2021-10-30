#ifndef __MCS_PIVOT_H__
#define __MCS_PIVOT_H__
#include "MCStep.h"

/*
Perhaps exclude the possibility to rotate te first triad.

*/

class MCS_Pivot;

class MCS_Pivot: public MCStep
{
protected:
    double fac_theta;
    double fac_size;
    arma::colvec  sigmas;

    std::uniform_int_distribution<> genhinge{0, 1};

public:
    MCS_Pivot(Chain * ch,const std::vector<long long int> & seedseq);
    MCS_Pivot(Chain * ch,const std::vector<long long int> & seedseq, int range_id1, int range_id2);
    ~MCS_Pivot();

    void update_settings();

// Monte Carlo Step (call main functionality)
public:
    bool MC_move();

};

#endif
