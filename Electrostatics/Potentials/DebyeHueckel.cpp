#include "DebyeHueckel.h"

DebyeHueckel::DebyeHueckel(Chain * ch, std::vector<double> & params, double integral_dx, double rho_max)
: ESPotential(ch,params,integral_dx,rho_max)
{
    if (params.size() != 2) {
        throw std::invalid_argument("DebyeHueckel::DebyeHueckel(): Debye Hückel potential requires 2 parameters.");
    }

    type                = "Debye Hückel";
    debye_length        = params[0];
    dielectric_constant = params[1];

//    prefactor = set_prefactor();
    prefactor = dielectric_constant;
}

DebyeHueckel::~DebyeHueckel() {
}

double DebyeHueckel::set_prefactor() {
    return 1;
}
double DebyeHueckel::integrant(double r) {
    /*
        The prefactor should not be set here. It is multiplied
        after summation.
    */
    return std::exp(-r/debye_length)/r;
}

