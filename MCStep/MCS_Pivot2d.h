#ifndef __MCS_PIVOT2D_H__
#define __MCS_PIVOT2D_H__
#include "MCStep.h"

/*
Perhaps exclude the possibility to rotate te first triad.

*/

class MCS_Pivot2d;

class MCS_Pivot2d: public MCStep
{
protected:
    double fac;
    double sigma;
    arma::colvec  normal;

    std::uniform_int_distribution<> genhinge{0, 1};

public:
    MCS_Pivot2d(Chain * ch,const std::vector<long long int> & seedseq, arma::colvec normal_to_plane);
    MCS_Pivot2d(Chain * ch,const std::vector<long long int> & seedseq, arma::colvec normal_to_plane, int range_id1, int range_id2);
    ~MCS_Pivot2d();

    void update_settings();

// Monte Carlo Step (call main functionality)
public:
    bool MC_move();


};

#endif
