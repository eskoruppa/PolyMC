#include "../BPStep.h"
#ifndef BPS_USE_EVALENERGY

std::vector<double> BPStep::stiffmat_set_params(std::vector<double> params) {
    double factor = 0.5/disc_len;
    for (unsigned i=0;i<params.size();i++) {
        params[i]*=REF_TEMP/temp*factor;
    }
    return params;
}

double BPStep::stiffmat_cal_beta_energy_diag(const arma::colvec & Theta) {
    return arma::dot(Theta, params_diag_mat * Theta);
}

double BPStep::stiffmat_cal_beta_energy_offdiag(const arma::colvec & Theta1,const arma::colvec & Theta2, const arma::mat M) {
    return arma::dot(Theta1, M * Theta2);
}

std::vector<double> BPStep::stiffmat_set_params_temp(double new_temp,std::vector<double> params) {
    for (unsigned i=0;i<params.size();i++) {
        params[i]*=temp/new_temp;
    }
    return params;
}



#endif


