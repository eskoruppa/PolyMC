#ifndef __MCS_STATESWITCH_H__
#define __MCS_STATESWITCH_H__
#include "MCStep.h"


/*
Perhaps exclude the possibility to rotate te first triad.

*/

class MCS_Stateswitch;

class MCS_Stateswitch: public MCStep
{
protected:
    bool   closed;
    int    rge_id1,rge_id2,rge_span;
    bool   closed_topol_restricted_range = false;
    arma::colvec interface_energies;
    std::uniform_int_distribution<> genhinge1{0, 1};
    std::uniform_int_distribution<> genhingedist{0, 1};

public:
    MCS_Stateswitch( Chain * ch,const std::vector<long long int> & seedseq, int selrange_min, int selrange_max);
    MCS_Stateswitch( Chain * ch,const std::vector<long long int> & seedseq, int selrange_min, int selrange_max, int range_id1, int range_id2);
    ~MCS_Stateswitch();

    void update_settings();
    void init_interface_energies();

// Monte Carlo Step (call main functionality)
public:
    bool MC_move();

};

#endif
