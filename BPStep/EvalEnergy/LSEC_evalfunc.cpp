#include "EnergyFuncs.h"

double cal_LSEC(double theta, double A, double theta0, int k) {
    if (theta < theta0) {
        return 0.5*A*theta*theta;
    }
    else {
        double lin = 0.5*A*theta*theta;
        double sub = 0.5*A*theta0*theta0 - 0.5*A/k * theta0*(theta0-M_PI)*(1 - std::pow( (theta-M_PI)/(theta0-M_PI) ,2*k) );
        if (0.5*A*theta0*theta0>sub) {
            cout << "something is wrong!" << endl;
            cout << "theta = " << theta << endl;
            cout << "lin = " << lin << endl;
            cout << "sub = " << sub << endl;
        }

        return 0.5*A*theta0*theta0 - 0.5*A/k * theta0*(theta0-M_PI)*(1 - std::pow( (theta-M_PI)/(theta0-M_PI) ,2*k) );
    }
}
