#include "BPStep.h"
#ifndef BPS_USE_EVALENERGY

#define BPS_METHOD_STIFFMAT_NAMEID  "stiffmat"
#define BPS_METHOD_STIFFMAT_ID      1

//#define STIFFONLY
unsigned BPStep::get_method_id(const std::string & method) {
    if (method==BPS_METHOD_STIFFMAT_NAMEID) {
        return BPS_METHOD_STIFFMAT_ID;
    }

    /*
        New Methods here
    */



    std::cout << std::endl << "Invalid interaction type '" + method + "'" << std::endl;
    std::cout << "Valid types:" << std::endl;
    std::cout << " - '" << BPS_METHOD_STIFFMAT_NAMEID << "'" << std::endl;
    /*
        List new methods
    */

    std::cout << std::endl << "Simulation terminated." << std::endl;
    std::exit(0);
}

std::vector<double> BPStep::set_params(unsigned method, std::vector<double> params) {

    if (method==BPS_METHOD_STIFFMAT_ID) {
        return stiffmat_set_params(params);
    }

    /*
        New Methods here
    */


}

bool BPStep::allow_params_averaging(unsigned method) {

    if (method==BPS_METHOD_STIFFMAT_ID) {
        return true;
    }

    /*
        New Methods here
    */


}



double BPStep::cal_beta_energy_diag(const arma::colvec & Theta) {
    #ifdef STIFFONLY
    stiffmat_cal_beta_energy_diag(Theta);
    #else

    if (coup_diag_method==BPS_METHOD_STIFFMAT_ID) {
        return stiffmat_cal_beta_energy_diag(Theta);
    }

    /*
        New Methods here
    */





    #endif
}


double BPStep::cal_beta_energy_left(const arma::colvec & Theta1,const arma::colvec & Theta2, unsigned id_left) {
    #ifdef STIFFONLY
    stiffmat_cal_beta_energy_left(Theta1,Theta2,id_left);
    #else

    if (coup_left_method[id_left]==BPS_METHOD_STIFFMAT_ID) {
        return stiffmat_cal_beta_energy_offdiag(Theta1,Theta2,params_left_mat.slice(id_left));
    }

    /*
        New Methods here
    */



    #endif
}

double BPStep::cal_beta_energy_right(const arma::colvec & Theta1,const arma::colvec & Theta2, unsigned id_right) {
    #ifdef STIFFONLY
    stiffmat_cal_beta_energy_left(Theta1,Theta2,id_left);
    #else

    if (coup_right_method[id_right]==BPS_METHOD_STIFFMAT_ID) {
        return stiffmat_cal_beta_energy_offdiag(Theta1,Theta2,params_right_mat.slice(id_right));
    }

    /*
        New Methods here
    */



    #endif
}


std::vector<double> BPStep::set_params_temp(double new_temp,unsigned method, std::vector<double> params) {

    if (method==BPS_METHOD_STIFFMAT_ID) {
        return stiffmat_set_params_temp(new_temp, params);
    }

    /*
        New Methods here
    */


}

bool BPStep::check_isotropic_bending() {

    if (coup_diag_method==BPS_METHOD_STIFFMAT_ID) {
        return ( params_diag_mat(0,0) == params_diag_mat(1,1) && params_diag_mat(0,1) == 0 && params_diag_mat(0,2) == 0 && params_diag_mat(1,2) == 0 );
    }

    /*
        New Methods here. Defining this method is not necessary. If it is not defined isotropy is determined by Monte Carlo.
    */



    return cal_MC_isotropic_bending();
}


arma::mat BPStep::cal_cov() {
    if (coup_diag_method==BPS_METHOD_STIFFMAT_ID) {
        return arma::inv(params_diag_mat/STIFFMAT_FACTOR);
    }


    /*
        New Methods here. Defining this method is not necessary. If it is not defined the covariance matrix is calculated by means of Monte Carlo.
    */
    return cal_MC_cov();
}




#endif
