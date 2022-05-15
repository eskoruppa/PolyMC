#ifndef __EXCLUDED_ESPOTENTIAL_DEBYEHUECKEL__
#define __EXCLUDED_ESPOTENTIAL_DEBYEHUECKEL__

#include "../ESPotential.h"

class DebyeHueckel;

class DebyeHueckel : public ESPotential
{

protected:

    double debye_length;
    double dielectric_constant;

public:

    DebyeHueckel(Chain * ch, std::vector<double> & params, double integral_dx, double rho_max);
    ~DebyeHueckel();

public:
    virtual double set_prefactor();
    virtual double integrant(double r);



};


#endif
