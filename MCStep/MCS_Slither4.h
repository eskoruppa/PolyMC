#ifndef __MCS_SLITHER4_H__
#define __MCS_SLITHER4_H__
#include "MCStep.h"


/*
Perhaps exclude the possibility to rotate te first triad.

*/

class MCS_Slither4;

class MCS_Slither4: public MCStep
{
protected:
    double fac_phi,fac_theta;
    double        sigma_phi;
    arma::colvec  sigma_theta;

    bool closed;
    std::uniform_int_distribution<> genhingeA{0, 1};
    std::uniform_int_distribution<> genhingedist{0, 1};

public:
    MCS_Slither4( Chain * ch,const std::vector<long long int> & seedseq, double seg_size_min, double seg_size_max);
    ~MCS_Slither4();

    void update_settings();

// Monte Carlo Step (call main functionality)
public:
    bool MC_move();
    bool rotation(int idA, int idB);


};

#endif
