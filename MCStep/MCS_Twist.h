#ifndef __MCS_TWIST_H__
#define __MCS_TWIST_H__
#include "MCStep.h"

/*
Perhaps exclude the possibility to rotate te first triad.

*/

class MCS_Twist;

class MCS_Twist: public MCStep
{
protected:
    double  sigfac;
    double  sigma;

    std::uniform_int_distribution<> genhinge{0, 1};

public:
    MCS_Twist(Chain * ch,const std::vector<long long int> & seedseq);
    MCS_Twist(Chain * ch,const std::vector<long long int> & seedseq, int range_id1, int range_id2);
    ~MCS_Twist();

    void update_settings();

// Monte Carlo Step (call main functionality)
public:
    bool MC_move();


};

#endif
