#ifndef __EE_KINKXY_INCLUDED__
#define __EE_KINKXY_INCLUDED__

#include "../EvalEnergy.h"

#define KINKXY_A1 0
#define KINKXY_Ay 1
#define KINKXY_Ax 2
#define KINKXY_Tk 3
#define KINKXY_Al 4
#define KINKXY_xl 5
#define KINKXY_hl 6
#define KINKXY_Ar 7
#define KINKXY_xr 8
#define KINKXY_hr 9

#define KINKXY_xcl          10
#define KINKXY_xcr          11
#define KINKXY_sinTk        12
#define KINKXY_cosTk        13
#define KINKXY_left_active  14
#define KINKXY_right_active 15


class EE_KinkXY;

class EE_KinkXY: public EvalEnergy
{
//////////////////////////////////////////////////////////////////////////////////
// member variables //////////////////////////////////////////////////////////////
protected:
    double A1,Ay,Ax,Tk,Al,xl,hl,Ar,xr,hr;
    double xcl,xcr; // parabola bundary
    double sinTk,cosTk;
    bool   left_active  = false;
    bool   right_active = false;

public:
////////////////////////////////////////////////////////////////////////////////////
/////// constructor / destructor ///////////////////////////////////////////////////
    EE_KinkXY(const std::vector<double> & params,double disc_len, double temp, bool is_diag);
    ~EE_KinkXY();

////////////////////////////////////////////////////////////////////////////////////
/////////////////   main   /////////////////////////////////////////////////////////
    double cal_beta_energy_diag(const arma::colvec & Theta);
    double cal_beta_energy_offdiag(const arma::colvec & Theta1, const arma::colvec & Theta2);

////////////////////////////////////////////////////////////////////////////////////
/////////////////  setters /////////////////////////////////////////////////////////
    void set_temp(double temp);
    void set_params(const std::vector<double> & new_params);

////////////////////////////////////////////////////////////////////////////////////
/////////////////  getters /////////////////////////////////////////////////////////
//    arma::mat get_cov();
    bool isotropic_bending();

    std::vector<double> get_status_diag(const arma::colvec & Theta);

//    bool has_left_kink();
//    bool has_right_kink();
//    int get_kink_state(const arma::colvec & Theta);

    arma::colvec get_theta_swap_left_kink(const arma::colvec & Theta);
    arma::colvec get_theta_swap_right_kink(const arma::colvec & Theta);


};

#endif
