#ifndef __MCS_TOPOLPERM_H__
#define __MCS_TOPOLPERM_H__
#include "MCStep.h"


class MCS_TopolPerm;

class MCS_TopolPerm: public MCStep
{
protected:
    bool closed;
    int trigger_every;
    int nbpm1;
    std::uniform_int_distribution<> direction_selection{0, 1};

public:
    MCS_TopolPerm( Chain * ch,const std::vector<long long int> & seedseq, int trigger_every = 1);
    ~MCS_TopolPerm();
    void update_settings();

// Monte Carlo Step (call main functionality)
public:
    bool MC_move();
    // bool rotation(int idA, int idB);

};

#endif
