#ifndef __MCS_PIVOTKINK_H__
#define __MCS_PIVOTKINK_H__
#include "MCStep.h"

/*
Perhaps exclude the possibility to rotate te first triad.

*/

class MCS_PivotKink;

class MCS_PivotKink: public MCStep
{
protected:
    double fac_theta;
    double fac_size;

    arma::colvec  sigmas;

    double A = 50;
    double D = 100;
    double B = 100;
    double x0 = 1;
    double h  = 4;
    double theta = 0;
    double xcon = 0.5;
    double costheta = 1;
    double sintheta = 0;

    arma::mat R0,R0_T;


    std::uniform_int_distribution<> genhinge{0, 1};

public:
    MCS_PivotKink(Chain * ch, double A, double D, double B, double x0, double h, double theta, arma::vec T0);
    ~MCS_PivotKink();

    void update_settings();

// Monte Carlo Step (call main functionality)
public:
    bool MC_move();
    double cal_energy(arma::colvec theta);
};

#endif
