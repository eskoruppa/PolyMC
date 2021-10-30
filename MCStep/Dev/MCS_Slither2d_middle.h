#ifndef __MCS_SLITHER2D_H__
#define __MCS_SLITHER2D_H__
#include "MCStep.h"


class MCS_Slither2d;

class MCS_Slither2d: public MCStep
{
protected:
    double fac;
    double sigma;
    arma::colvec  normal;
    double max_rot;

    bool   closed;
    std::uniform_int_distribution<> genhingeA{0, 1};
    std::uniform_int_distribution<> genhingedist{0, 1};

public:
    MCS_Slither2d( Chain * ch, double seg_size_min, double seg_size_max, arma::colvec normal_to_plane);
    ~MCS_Slither2d();

    void update_settings();

// Monte Carlo Step (call main functionality)
public:
    bool MC_move();
    bool rotation(int idA, int idB);

    arma::colvec max_angle(arma::colvec& v1, arma::colvec& w, double lv23, double max_select);


};

#endif
