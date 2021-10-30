#ifndef __MCS_CSROT2D_H__
#define __MCS_CSROT2D_H__
#include "MCStep.h"


/*
Perhaps exclude the possibility to rotate te first triad.

*/

class MCS_CSrot2d;

class MCS_CSrot2d: public MCStep
{
protected:
    arma::colvec normal;
    bool         closed;

    std::uniform_int_distribution<> genhinge1{0, 1};
    std::uniform_int_distribution<> genhingedist{0, 1};

public:
    MCS_CSrot2d( Chain * ch,const std::vector<long long int> & seedseq, int selrange_min, int selrange_max,arma::colvec normal = {0,0,1});
    ~MCS_CSrot2d();

    void update_settings();

// Monte Carlo Step (call main functionality)
public:
    bool MC_move();
};

#endif
