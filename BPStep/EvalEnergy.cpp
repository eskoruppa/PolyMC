#include "EvalEnergy.h"

////////////////////////////////////////////////////////////////////////////////////
/////// constructor / destructor ///////////////////////////////////////////////////

EvalEnergy::EvalEnergy(const std::vector<double> & params,double disc_len, double temp, bool is_diag)
: params(params), disc_len(disc_len),temp(temp), is_diag(is_diag) {
    kBT = REF_KBT*temp/REF_TEMP;
    ref_temp = REF_TEMP;

}
EvalEnergy::~EvalEnergy() {}

////////////////////////////////////////////////////////////////////////////////
///////////////   main   ///////////////////////////////////////////////////////

double EvalEnergy::cal_beta_energy_diag(const arma::colvec & Theta) {
    throw std::logic_error("cal_beta_energy_diag() not defined in class inheriting EvalEnergy");
    return 0;
}

double EvalEnergy::cal_beta_energy_offdiag(const arma::colvec & Theta1, const arma::colvec & Theta2) {
    throw std::logic_error("cal_beta_energy_offdiag() not defined in class inheriting EvalEnergy");
    return 0;
}

////////////////////////////////////////////////////////////////////////////////
///////////////  setters ///////////////////////////////////////////////////////
void EvalEnergy::set_temp(double temp) {
    this->temp = temp;
    this->kBT  = REF_KBT*temp/REF_TEMP;
}

void EvalEnergy::set_params(const std::vector<double> & params) {
    throw std::logic_error("set_params() not defined in class inheriting EvalEnergy");
}

void EvalEnergy::deactivate_twist_energy() {
    twist_energy_active = false;
}
void EvalEnergy::reactivate_twist_energy() {
    twist_energy_active = true;
}

////////////////////////////////////////////////////////////////////////////////
///////////////  getters ///////////////////////////////////////////////////////
arma::mat EvalEnergy::get_cov() {
    return cal_MC_cov();
}
std::string * EvalEnergy::get_method() {
    return &method;
}
std::vector<double> EvalEnergy::get_params() {
    return params;
}
std::vector<double> * EvalEnergy::get_attributes(){
    return &attributes;
}

bool EvalEnergy::get_parameter_averaging_allowed() {
    return parameter_averaging_allowed;
}

bool EvalEnergy::isotropic_bending() {

    std::cout << "isotropic_bending starting" << std::endl;

    std::random_device                      rd{};
    std::mt19937                            gen{rd()};
    std::uniform_real_distribution<double>  uniform{0.0,2*M_PI};
    std::seed_seq seed{rd(), rd(), rd(), rd(), rd(), rd(), rd(), rd()};
    gen.seed(seed);

    double E1,E2;
    double phi,cosphi,sinphi;
    arma::colvec Theta1,Theta2;
    Theta2 = arma::zeros(3);
    for (long long int i=0;i<ISOTROPIC_BENDING_NUM_SAMPLED;i++){
        Theta1    = {uniform(gen),uniform(gen),uniform(gen)};
        phi      = uniform(gen);
        sinphi   = std::sin(phi);
        cosphi   = std::cos(phi);
        E1       = cal_beta_energy_diag(Theta1);
        Theta2(0) = Theta1(0)*cosphi-Theta1(1)*sinphi;
        Theta2(1) = Theta1(0)*sinphi+Theta1(1)*cosphi;
        Theta2(2) = Theta1(2);
        E2       = cal_beta_energy_diag(Theta2);
        if (std::abs(E2-E1)>ISOTROPIC_BENDING_EQUAL_THRES) {
            std::cout << std::abs(E2-E1) << std::endl;
            std::cout << "Is NOT isotropic!" << std::endl;
            return false;
        }
    }
    std::cout << "Is isotropic!" << std::endl;
    return true;
}

std::vector<double> EvalEnergy::get_status_diag(const arma::colvec & Theta) {
    return {};
}

////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

arma::mat EvalEnergy::cal_MC_cov(long long int steps,double sigma) {
    std::random_device                      rd{};
    std::mt19937                            gen{rd()};
    std::normal_distribution<double>        normaldist{0.0,1.0};
    std::uniform_real_distribution<double>  uniformdist{0.0,1.0};
    std::seed_seq seed{rd(), rd(), rd(), rd(), rd(), rd(), rd(), rd()};
    gen.seed(seed);

    arma::mat    cov    = arma::zeros(3,3);
    arma::colvec Theta    = arma::zeros(3);
    arma::colvec newTheta = arma::zeros(3);
    double currentbetaE,newbetaE;
    currentbetaE = cal_beta_energy_diag(Theta);
    for (long long int i=0;i<steps;i++){
        newTheta(0) = Theta(0) + normaldist(gen)*sigma;
        newTheta(1) = Theta(1) + normaldist(gen)*sigma;
        newTheta(2) = Theta(2) + normaldist(gen)*sigma;
        newbetaE = cal_beta_energy_diag(newTheta);
        if (std::exp( currentbetaE - newbetaE ) > uniformdist(gen)) {
            Theta         = newTheta;
            currentbetaE  = newbetaE;
        }
        cov    += Theta*Theta.t();
    }
    return cov /= steps;
}
