//// include input functions
#include "Input/Argparse.h"
#include "Input/InputRead.h"

#include "PolyMC.h"
#include "MCMethods/ReplicaExchange.h"
#define POLYMC_VERSION 0.610

/*
TO:
    - Tweezer Torque requires sigma to be set to zero upon initialization for dLK to stay consistent (not sure why).



*/


// Defines
/*
ARMA_DONT_USE_WRAPPER
ARMA_USE_LAPACK
ARMA_USE_BLAS
ARMA_USE_CXX11

ARMA_NO_DEBUG
*/

// linkers
/*
-llapack

-lblas
*/

// Other Linker Options
/*
-larmadillo
*/

#include "Input/ManipulateArgv.h"


int main(int argc, const char **argcv) {
    std::cout << "PolyMC version: " << POLYMC_VERSION << std::endl;
    arma::arma_version ver;
    std::cout << "ARMA version: "<< ver.as_string() << std::endl;

    std::vector<std::string> argv(argcv, argcv+argc);

    std::string inputfn = "";
    inputfn  = parse_arg(inputfn, "-in",argv);
    InputRead * input = new InputRead(inputfn);

    if (input->contains_multi("replica_exchange")) {
        ReplicaExchange re(inputfn,argv);
    }
    else {
        PolyMC polyMC(argv);
        polyMC.simple_execute();
    }

    return 0;
}
