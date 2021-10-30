#ifndef __EVALENERGY_INCLUDED__
#define __EVALENERGY_INCLUDED__

#include <armadillo> // arma
#include <iostream> // std::cout
#include <string>   // std::string
#include <algorithm>
#include <stdexcept>
#include <vector>
#include <cmath>
#include <random>

#include "../ExtraFuncs.h"
#include "../SO3Methods.h"


#define REF_TEMP 300
#define REF_KBT  4.114

#define MC_COV_DEFAULT_STEPS 1e6
//#define MC_COV_DEFAULT_STEPS 1e4
#define MC_COV_DEFAULT_SIGMA 0.1

#define ISOTROPIC_BENDING_NUM_SAMPLED 2e5
#define ISOTROPIC_BENDING_EQUAL_THRES 1e-8


class EvalEnergy;

class EvalEnergy {

//////////////////////////////////////////////////////////////////////////////////
// member variables //////////////////////////////////////////////////////////////
protected:
    std::string         method="template";
    bool parameter_averaging_allowed = false;

    std::vector<double> params;
    std::vector<double> attributes;
    double disc_len;    // discretization length
    double temp;        // temperature
    double ref_temp;    // Should be 300 Kelvin. Stiffness Matrix should be expressed in terms of units of beta (1/(4.114))
    double kBT;         // 4.114pNnm by default (at 300K)
    bool   is_diag;     // couples the same Theta (true) or two different Thetas

    double twist_energy_active = true;

public:
//////////////////////////////////////////////////////////////////////////////////
/////// constructor / destructor /////////////////////////////////////////////////
    EvalEnergy(const std::vector<double> & params,double disc_len, double temp, bool is_diag);
    virtual ~EvalEnergy();

////////////////////////////////////////////////////////////////////////////////
///////////////   main   ///////////////////////////////////////////////////////
    virtual double cal_beta_energy_diag(const arma::colvec & Theta);
    virtual double cal_beta_energy_offdiag(const arma::colvec & Theta1, const arma::colvec & Theta2);

////////////////////////////////////////////////////////////////////////////////
///////////////  setters ///////////////////////////////////////////////////////
    virtual void set_temp(double temp);
    virtual void set_params(const std::vector<double> & new_params);

    virtual void deactivate_twist_energy();
    virtual void reactivate_twist_energy();

////////////////////////////////////////////////////////////////////////////////
///////////////  getters ///////////////////////////////////////////////////////
    virtual arma::mat get_cov();

    std::string * get_method();
    std::vector<double> get_params();
    std::vector<double> * get_attributes();

    bool get_parameter_averaging_allowed();

    virtual bool isotropic_bending();

    virtual std::vector<double> get_status_diag(const arma::colvec & Theta);

////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
protected:
    arma::mat cal_MC_cov(long long int steps=MC_COV_DEFAULT_STEPS, double sigma=MC_COV_DEFAULT_SIGMA);

};

#endif
