#ifndef __MCS_CSROT_H__
#define __MCS_CSROT_H__
#include "MCStep.h"


/*
Perhaps exclude the possibility to rotate te first triad.

*/

class MCS_CSrot;

class MCS_CSrot: public MCStep
{
protected:
    double fac;
    double sigma;
    bool   closed;
    int    rge_id1,rge_id2,rge_span;
    bool   closed_topol_restricted_range = false;
    std::uniform_int_distribution<> genhinge1{0, 1};
    std::uniform_int_distribution<> genhingedist{0, 1};

public:
    MCS_CSrot( Chain * ch,const std::vector<long long int> & seedseq, int selrange_min, int selrange_max);
    MCS_CSrot( Chain * ch,const std::vector<long long int> & seedseq, int selrange_min, int selrange_max, int range_id1, int range_id2);
    ~MCS_CSrot();

    void update_settings();

// Monte Carlo Step (call main functionality)
public:
    bool MC_move();
    bool rotation(int idA, int idB);

protected:
    virtual double gen_trial_conf() override;


};

#endif
