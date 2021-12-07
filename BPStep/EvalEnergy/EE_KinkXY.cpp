#include "EE_KinkXY.h"

////////////////////////////////////////////////////////////////////////////////////
/////// constructor / destructor ///////////////////////////////////////////////////
EE_KinkXY::EE_KinkXY(const std::vector<double> & params, double disc_len, double temp, bool is_diag)
: EvalEnergy(params, disc_len, temp, is_diag)
{
    if (!is_diag) {
        std::cout << "Error: KinkXY can only be used for on-site couplings!" << std::endl;
        std::exit(0);
    }

    method = "kinkxy";
    parameter_averaging_allowed = false;
    set_params(params);
}

EE_KinkXY::~EE_KinkXY() {}

////////////////////////////////////////////////////////////////////////////////////
/////////////////   main   /////////////////////////////////////////////////////////
double EE_KinkXY::cal_beta_energy_diag(const arma::colvec & Theta) {
    double E,x,y;
    x =  cosTk* Theta(1) - sinTk* Theta(2);
    y =  sinTk* Theta(1) + cosTk* Theta(2);

    // Theta1 contribution
    E = 0.5*A1*Theta(0)*Theta(0);

    // y contribution
    E += 0.5*Ay*y*y;

    // x contribution
    if (left_active && x<xcl) {
        E += 0.5*Al*(x-xl)*(x-xl) + hl;
        return E;
    }
    if (right_active && x>xcr) {
        E += 0.5*Ar*(x-xr)*(x-xr) + hr;
        return E;
    }
    E += 0.5*Ax*x*x;
    return E;
}

double EE_KinkXY::cal_beta_energy_offdiag(const arma::colvec & Theta1, const arma::colvec & Theta2) {
    std::cout << "Error: KinkXY can only be used for on-site couplings!" << std::endl;
    std::exit(0);
    return 0;
}

////////////////////////////////////////////////////////////////////////////////////
/////////////////  setters /////////////////////////////////////////////////////////
void EE_KinkXY::set_temp(double temp) {
    A1 = A1 * this->temp/temp;
    Ay = Ay * this->temp/temp;
    Ax = Ax * this->temp/temp;
    Al = Al * this->temp/temp;
    hl = hl * this->temp/temp;
    Ar = Ar * this->temp/temp;
    hr = hr * this->temp/temp;
    EvalEnergy::set_temp(temp);

    attributes[KINKXY_A1] = A1;
    attributes[KINKXY_Ay] = Ay;
    attributes[KINKXY_Ax] = Ax;
    attributes[KINKXY_Al] = Al;
    attributes[KINKXY_hl] = hl;
    attributes[KINKXY_Ar] = Ar;
    attributes[KINKXY_hr] = hr;
}


void EE_KinkXY::set_params(const std::vector<double> & params) {
    if ( params.size() < 10 ) {
        throw std::invalid_argument("Insufficient parameters provided for kinkxy (at least 10 required)!");
    }
    this->params = params;

    A1 = params[KINKXY_A1] / disc_len * REF_TEMP/temp;
    Ay = params[KINKXY_Ay] / disc_len * REF_TEMP/temp;
    Ax = params[KINKXY_Ax] / disc_len * REF_TEMP/temp;
    Tk = params[KINKXY_Tk] / 180 * M_PI;
    Al = params[KINKXY_Al] / disc_len * REF_TEMP/temp;
    xl = params[KINKXY_xl] / 180 * M_PI;
    hl = params[KINKXY_hl] * REF_TEMP/temp;
    Ar = params[KINKXY_Ar] / disc_len * REF_TEMP/temp;
    xr = params[KINKXY_xr] / 180 * M_PI;
    hr = params[KINKXY_hr] * REF_TEMP/temp;

    sinTk = std::sin(Tk);
    cosTk = std::cos(Tk);

/*
 *         Kl             |
 *          \             |
 *           \     _      |          Kr
 *            \   / \     |     _   /
 *   hr .....  \_/ | \    |    / \_/ ..... hr
 *             xl  |  \ Kx|   / | xr
 *                 |   \  |  /  |
 * ----------------|------|-----|------------------
 *                xcl          xcr
 */

    if (A1<0) {
        throw std::invalid_argument("Invlid input for KinkXY: A1 < 0!");
    }
    if (Ay<0) {
        throw std::invalid_argument("Invlid input for KinkXY: Ay < 0!");
    }
    if (Ax<0) {
        throw std::invalid_argument("Invlid input for KinkXY: Ax < 0!");
    }
    if (Al<0) {
        throw std::invalid_argument("Invlid input for KinkXY: Al < 0!");
    }
    if (Ar<0) {
        throw std::invalid_argument("Invlid input for KinkXY: Ar < 0!");
    }
    if (xl>0) {
        throw std::invalid_argument("Invlid input for KinkXY: xl > 0!");
    }
    if (xr<0) {
        throw std::invalid_argument("Invlid input for KinkXY: xr < 0!");
    }

    left_active = false;
    if (xl < 0) {
        double discr = Ax*Al*xl*xl + 2*hl*(Ax-Al);
        if (discr > 0) {
            left_active = true;
            if (Ax == Al) {
                xcl = 0.5*xl+hl/(Al*xl);
            }
            else{
                xcl = (-Al*xl + std::sqrt(discr)) / (Ax - Al);
                if (xcl > 0 ) {
                    xcl = (-Al*xl - std::sqrt(discr)) / (Ax - Al);
                }
                if (xcl > 0 ) {
                    left_active = false;
                }
            }
        }
    }

    right_active = false;
    if (xr > 0) {
        double discr = Ax*Ar*xr*xr + 2*hr*(Ax-Ar);
        if (discr > 0) {
            right_active = true;
            if (Ax == Ar) {
                xcr = 0.5*xr+hr/(Ar*xr);
            }
            else{
                xcr = (-Ar*xr - std::sqrt(discr)) / (Ax - Ar);
                if (xcr < 0) {
                    xcr = (-Ar*xr + std::sqrt(discr)) / (Ax - Ar);
                }
                if (xcr < 0) {
                    right_active = false;
                }
            }
        }
    }

    // set attributes
    attributes.clear();
    attributes.push_back(A1);
    attributes.push_back(Ay);
    attributes.push_back(Ax);
    attributes.push_back(Tk);
    attributes.push_back(Al);
    attributes.push_back(xl);
    attributes.push_back(hl);
    attributes.push_back(Ar);
    attributes.push_back(xr);
    attributes.push_back(hr);
    attributes.push_back(xcl);
    attributes.push_back(xcr);
    attributes.push_back(sinTk);
    attributes.push_back(cosTk);
    attributes.push_back(left_active);
    attributes.push_back(right_active);
}



////////////////////////////////////////////////////////////////////////////////////
/////////////////  getters /////////////////////////////////////////////////////////

bool EE_KinkXY::isotropic_bending() {
    if (left_active || right_active) {
        return false;
    }
    if (Tk==0 && Ax == A1) {
        return true;
    }
    return false;
}

std::vector<double> EE_KinkXY::get_status_diag(const arma::colvec & Theta) {
/*
    returns length 1 vector indicating whether bps is in kinked state
    unkinked:  0
    left:     -1
    right:    +1
*/
    double kinked_state = 0;

    double x;
    x =  cosTk* Theta(1) - sinTk* Theta(2);

    // x contribution
    if (left_active && x<xcl) {
        kinked_state = -1;
    }
    if (right_active && x>xcr) {
        kinked_state =  1;
    }
    return {kinked_state,double(left_active),double(right_active)};
}


//bool EE_KinkXY::has_left_kink() {
//    return left_active;
//}
//bool EE_KinkXY::has_right_kink() {
//    return right_active;
//}
//
//int EE_KinkXY::get_kink_state(const arma::colvec & Theta) {
///*
//    returns length 1 vector indicating whether bps is in kinked state
//    unkinked:  0
//    left:     -1
//    right:    +1
//*/
//    double x;
//    x =  cosTk* Theta(1) - sinTk* Theta(2);
//    // x contribution
//    if (left_active && x<xcl) {
//        return -1;
//    }
//    if (right_active && x>xcr) {
//        return 1;
//    }
//    return 0;
//}

arma::colvec EE_KinkXY::get_theta_swap_left_kink(const arma::colvec & Theta) {
    if (!left_active) {
        // if left kink is not active return current Theta
        return Theta;
    }

    double x,y;
    x =  cosTk* Theta(1) - sinTk* Theta(2);
    y =  sinTk* Theta(1) + cosTk* Theta(2);

    if (right_active && x>xcr) {
        // if BPS is currently in right kink return current Theta (no move)
        return Theta;
    }

    double x_new;
    if (x<xcl) {
        // is currently in left kink
        x_new = x - xcl;
    }
    else {
        // is currently not kinked
        x_new = x + xcl;
    }

    arma::colvec newTheta = arma::zeros(3);
    newTheta(1) =  cosTk* x_new + sinTk* y;
    newTheta(2) = -sinTk* x_new + cosTk* y;
    return newTheta;
}

arma::colvec EE_KinkXY::get_theta_swap_right_kink(const arma::colvec & Theta) {
    if (!right_active) {
        // if right kink is not active return current Theta
        return Theta;
    }

    double x,y;
    x =  cosTk* Theta(1) - sinTk* Theta(2);
    y =  sinTk* Theta(1) + cosTk* Theta(2);

    if (left_active && x<xcl) {
        // if BPS is currently in left kink return current Theta (no move)
        return Theta;
    }

    double x_new;
    if (x>xcr) {
        // is currently in right kink
        x_new = x - xcr;
    }
    else {
        // is currently not kinked
        x_new = x + xcr;
    }

    arma::colvec newTheta = arma::zeros(3);
    newTheta(1) =  cosTk* x_new + sinTk* y;
    newTheta(2) = -sinTk* x_new + cosTk* y;
    return newTheta;
}
