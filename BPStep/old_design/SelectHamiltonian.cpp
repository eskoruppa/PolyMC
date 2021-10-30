#include "BPStep.h"



EvalEnergy * BPStep::select_EvalEnergy(const std::string & method, const std::vector<double> & params, bool is_diag) {

    if (method=="stiffmat") {
        EE_StiffMat * stiffmat = new EE_StiffMat(params,disc_len,T,is_diag);
        return stiffmat;
    }


    /*
        Add new methods here
    */



    std::cout << std::endl << "Invalid interaction type '" + method + "'" << std::endl;
    std::cout << "Valid types:" << std::endl;
    std::cout << " - 'stiffmat'" << std::endl;

    /*
        List new methods
    */

    std::cout << std::endl << "Simulation terminated." << std::endl;
    std::exit(0);
}


double BPStep::eval(const arma::colvec & Theta) {
    return arma::dot(Theta, M0 * Theta);
}

