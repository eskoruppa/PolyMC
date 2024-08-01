#include "EE_Bubble.h"

////////////////////////////////////////////////////////////////////////////////////
/////// constructor / destructor ///////////////////////////////////////////////////

EE_Bubble::EE_Bubble(const std::vector<double> & params,double disc_len, double temp, bool is_diag)
: EvalEnergy(params, disc_len, temp, is_diag)
{
    if (!is_diag) {
        std::cout << "Error: Bubble can only be used for on-site couplings!" << std::endl;
        std::exit(0);
    }

    method = "bubble";
    parameter_averaging_allowed = false;

    set_params(params);
}

EE_Bubble::~EE_Bubble() {}


////////////////////////////////////////////////////////////////////////////////////
/////////////////   stateswitch   //////////////////////////////////////////////////

void EE_Bubble::propose_stateswitch(int new_state) {
    if (new_state < 1 || new_state > 2) {
        std::cout << "Error: invalid state " << new_state << " for EE_Bubble." << std::endl;
        std::exit(0); 
    }

    if (new_state == 1) {
        proposed_state = 1;
        proposed_stiffmat = state1_stiffmat;
        proposed_stateenergy = state1_stateenergy;
    }
    else {
        proposed_state = 2;
        proposed_stiffmat = state2_stiffmat;
        proposed_stateenergy = state2_stateenergy;
    }
    stateswitch_pending = true;
}

void EE_Bubble::set_switch(bool accepted) {
    if (accepted) {
        current_state = proposed_state;
        current_stiffmat = proposed_stiffmat;
        current_stateenergy = proposed_stateenergy;
    }
    stateswitch_pending = false;
}

////////////////////////////////////////////////////////////////////////////////////
/////////////////   main   /////////////////////////////////////////////////////////
double EE_Bubble::cal_beta_energy_diag(const arma::colvec & Theta) {
    if (stateswitch_pending) {
        return arma::dot(Theta, proposed_stiffmat * Theta) + proposed_stateenergy;
    }
    return arma::dot(Theta, current_stiffmat * Theta) + current_stateenergy;
}

double EE_Bubble::cal_beta_energy_offdiag(const arma::colvec & Theta1, const arma::colvec & Theta2) {
    std::cout << "Error: Bubble can only be used for on-site couplings!" << std::endl;
    std::exit(0);
}

////////////////////////////////////////////////////////////////////////////////////
/////////////////  setters /////////////////////////////////////////////////////////
void EE_Bubble::set_temp(double temp) {

    state1_stiffmat = state1_stiffmat *this->temp/temp;
    state2_stiffmat = state2_stiffmat *this->temp/temp;
    state1_stateenergy = state1_stateenergy *this->temp/temp;
    state2_stateenergy = state2_stateenergy *this->temp/temp;

    current_stiffmat = current_stiffmat *this->temp/temp;
    current_stateenergy = current_stateenergy *this->temp/temp;

    interface_energy = interface_energy *this->temp/temp;

    EvalEnergy::set_temp(temp);
}

void EE_Bubble::set_params(const std::vector<double> & params) {
    if ( params.size() < 21 ) {
        throw std::invalid_argument("Insufficient parameters provided for stiffness matrix (at least 20 required)!");
    }
    this->params = params;
    factor = 0.5;
    state1_stiffmat = {{params[0],params[1],params[2]},{params[3],params[4],params[5]},{params[6],params[7],params[8]}};
    state1_stiffmat *= factor/disc_len * REF_TEMP/temp;
    state1_stateenergy = params[9] * REF_TEMP/temp;

    state2_stiffmat = {{params[10],params[11],params[12]},{params[13],params[14],params[15]},{params[16],params[17],params[18]}};
    state2_stiffmat *= factor/disc_len * REF_TEMP/temp;
    state2_stateenergy = params[19] * REF_TEMP/temp;

    interface_energy = params[20] * REF_TEMP/temp;

    current_state = 1;
    stateswitch_pending = false;

    current_stiffmat = state1_stiffmat;
    current_stateenergy = state1_stateenergy;
}

void EE_Bubble::deactivate_twist_energy() {
    std::cout << "Error: Bubble does not allow for deactivation of twist energy!" << std::endl;
    std::exit(0);
}

void EE_Bubble::reactivate_twist_energy() {
    std::cout << "Error: Bubble does not allow for deactivation of twist energy!" << std::endl;
    std::exit(0);
}

////////////////////////////////////////////////////////////////////////////////////
/////////////////  getters /////////////////////////////////////////////////////////
arma::mat EE_Bubble::get_cov() {
    return arma::inv(state1_stiffmat/factor);
}

bool EE_Bubble::isotropic_bending() {
    return ( state1_stiffmat(0,0) == state1_stiffmat(1,1) && state1_stiffmat(0,1) == 0 && state1_stiffmat(0,2) == 0 && state1_stiffmat(1,2) == 0 );
}

double EE_Bubble::get_interface_energy() {
    return interface_energy;
}