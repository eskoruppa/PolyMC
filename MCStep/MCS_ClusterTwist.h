#ifndef __MCS_CLUSTERTWIST_H__
#define __MCS_CLUSTERTWIST_H__
#include "MCStep.h"


/*
Perhaps exclude the possibility to rotate te first triad.

*/

class MCS_ClusterTwist;

class MCS_ClusterTwist: public MCStep
{
protected:
    double sigfac;
    double sigma;
    int    rge_id1,rge_id2,rge_span;
    bool   closed_topol_restricted_range = false;

    bool quick_eval;
    bool closed;

    std::uniform_int_distribution<> genhinge1{0, 1};
    std::uniform_int_distribution<> genhingedist{0, 1};

public:
    MCS_ClusterTwist( Chain * ch,const std::vector<long long int> & seedseq, int selrange_min, int selrange_max);
    MCS_ClusterTwist( Chain * ch,const std::vector<long long int> & seedseq, int selrange_min, int selrange_max, int range_id1, int range_id2);
    ~MCS_ClusterTwist();

    void update_settings();

// Monte Carlo Step (call main functionality)
public:
    bool MC_move();
    bool rotation(int idA, int idB);


};

#endif
