#include "EE_StiffMat.h"

////////////////////////////////////////////////////////////////////////////////////
/////// constructor / destructor ///////////////////////////////////////////////////

EE_StiffMat::EE_StiffMat(const std::vector<double> & params,double disc_len, double temp, bool is_diag)
: EvalEnergy(params, disc_len, temp, is_diag)
{
    method = "stiffmat";
    parameter_averaging_allowed = true;

    set_params(params);
}

EE_StiffMat::~EE_StiffMat() {}

////////////////////////////////////////////////////////////////////////////////////
/////////////////   main   /////////////////////////////////////////////////////////
double EE_StiffMat::cal_beta_energy_diag(const arma::colvec & Theta) {
    return arma::dot(Theta, current_stiffmat * Theta);
}

double EE_StiffMat::cal_beta_energy_offdiag(const arma::colvec & Theta1, const arma::colvec & Theta2) {
    return arma::dot(Theta1, current_stiffmat * Theta2);
}

////////////////////////////////////////////////////////////////////////////////////
/////////////////  setters /////////////////////////////////////////////////////////
void EE_StiffMat::set_temp(double temp) {
    stiffmat = stiffmat *this->temp/temp;
    set_current_stiffmat();
    EvalEnergy::set_temp(temp);
}

void EE_StiffMat::set_params(const std::vector<double> & params) {
    if ( params.size() < 9 ) {
        throw std::invalid_argument("Insufficient parameters provided for stiffness matrix (at least 9 required)!");
    }
    this->params = params;
    factor = 0.5;
//    if (!is_diag) factor = 0.25;
    stiffmat = {{params[0],params[1],params[2]},{params[3],params[4],params[5]},{params[6],params[7],params[8]}};
    stiffmat *= factor/disc_len * REF_TEMP/temp;
    set_current_stiffmat();
}

void EE_StiffMat::set_current_stiffmat() {
    current_stiffmat = stiffmat;
    if (!twist_energy_active) {
        current_stiffmat(2,2) = 0;
        current_stiffmat(2,1) = 0;
        current_stiffmat(1,2) = 0;
        current_stiffmat(2,0) = 0;
        current_stiffmat(0,2) = 0;
    }
}

void EE_StiffMat::deactivate_twist_energy() {
    twist_energy_active = false;
    set_current_stiffmat();
}
void EE_StiffMat::reactivate_twist_energy() {
    twist_energy_active = true;
    set_current_stiffmat();
}


////////////////////////////////////////////////////////////////////////////////////
/////////////////  getters /////////////////////////////////////////////////////////
arma::mat EE_StiffMat::get_cov() {
    return arma::inv(stiffmat/factor);
}

bool EE_StiffMat::isotropic_bending() {
    return ( stiffmat(0,0) == stiffmat(1,1) && stiffmat(0,1) == 0 && stiffmat(0,2) == 0 && stiffmat(1,2) == 0 );
}


// ////////////////////////////////////////////////////////////////////////////////
// /////////////////   stateswitch   //////////////////////////////////////////////
// void EE_StiffMat::propose_stateswitch(int new_state) {
//     throw std::logic_error("propose_stateswitch() not defined in class EE_KinkXY");
// }

// void EE_StiffMat::set_switch(bool accepted) {
//     throw std::logic_error("set_switch() not defined in class class EE_KinkXY");
// }