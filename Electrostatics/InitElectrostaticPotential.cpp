#include "../PolyMC.h"
#include "ESPotential.h"
#include "Potentials/DebyeHueckel.h"

bool PolyMC::init_electrostatic_potential(const std::string & type, std::vector<double> & params, double integral_dx, double rho_max){

    if (type == "debye_hueckel") {
        DebyeHueckel * dehu = new DebyeHueckel(chain,params,integral_dx, rho_max);
        ES_Potential = dehu;
        return true;
    }

    throw std::invalid_argument("init_electrostatic_potential(): Unknown electrostatics potential type");



    return false;
}




