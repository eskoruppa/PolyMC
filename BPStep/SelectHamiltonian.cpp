#include "BPStep.h"
#ifdef BPS_USE_EVALENERGY


EvalEnergy * BPStep::select_EvalEnergy(const std::string & method, const std::vector<double> & params, bool is_diag) {

    if (method=="stiffmat") {
        EE_StiffMat * stiffmat = new EE_StiffMat(params,disc_len,temp,is_diag);
        return stiffmat;
    }

    if (method=="kinkxy") {
        EE_KinkXY * kinkxy = new EE_KinkXY(params,disc_len,temp,is_diag);
        return kinkxy;
    }

    if (method=="bubble") {
        EE_Bubble * bubble = new EE_Bubble(params,disc_len,temp,is_diag);
        return bubble;
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

#endif
