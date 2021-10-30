#ifndef __MCS_KINKXY_H__
#define __MCS_KINKXY_H__
#include "MCStep.h"

/*
Perhaps exclude the possibility to rotate te first triad.

*/

class MCS_Kinkxy;

class MCS_Kinkxy: public MCStep
{
protected:
    double fac_sigma;
    double sigma;
    int num_in_tail;
//    arma::colvec  sigmas;
    std::uniform_int_distribution<> genhinge{0, 1};




public:
    MCS_Kinkxy(Chain * ch,const std::vector<long long int> & seedseq);
    ~MCS_Kinkxy();

    void update_settings();

// Monte Carlo Step (call main functionality)
public:
    bool MC_move();

protected:
    bool left_tail_move(int hingeID);
    bool right_tail_move(int hingeID);
    bool double_tail_move(int hingeID);


    bool get_kink_swap_Theta(int bpsid, const arma::colvec & Theta, arma::colvec & newTheta);

    int get_kink_state(const arma::colvec & Theta , const std::vector<double> & attr);

    bool get_theta_swap_left_kink  (const arma::colvec & Theta, const std::vector<double> & attr, arma::colvec & newTheta);
    bool get_theta_swap_right_kink (const arma::colvec & Theta, const std::vector<double> & attr, arma::colvec & newTheta);

    // test reversibility
    bool test_reverse(int hingeID);

};

#endif
