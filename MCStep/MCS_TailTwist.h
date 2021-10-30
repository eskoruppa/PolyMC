#ifndef __MCS_TAILTWIST_H__
#define __MCS_TAILTWIST_H__
#include "MCStep.h"

/*
Perhaps exclude the possibility to rotate te first triad.

*/

class MCS_TailTwist;

class MCS_TailTwist: public MCStep
{
protected:
    double sigfac;
    double sigma;

    bool quick_eval;
    std::uniform_int_distribution<> genhinge{0, 1};

public:
    MCS_TailTwist(Chain * ch,const std::vector<long long int> & seedseq);
    MCS_TailTwist(Chain * ch,const std::vector<long long int> & seedseq, int range_id1);
    ~MCS_TailTwist();

    void update_settings();

// Monte Carlo Step (call main functionality)
public:
    bool MC_move();


};

#endif
