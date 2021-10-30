#ifndef __MCS_ZPIV_H__
#define __MCS_ZPIV_H__
#include "MCStep.h"

/*
Perhaps exclude the possibility to rotate te first triad.

*/

class MCS_zPiv;

class MCS_zPiv: public MCStep
{
protected:
    double fac_sigma;
    double sigma;
//    arma::colvec  sigmas;
    std::uniform_int_distribution<> genhinge{0, 1};

public:
    MCS_zPiv(Chain * ch,const std::vector<long long int> & seedseq);
    MCS_zPiv(Chain * ch,const std::vector<long long int> & seedseq, int range_id1, int range_id2);
    ~MCS_zPiv();

    void update_settings();

// Monte Carlo Step (call main functionality)
public:
    bool MC_move();


};

#endif
