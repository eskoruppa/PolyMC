//// include input functions
#include "Input/Argparse.h"
#include "Input/InputRead.h"
#include "Input/ManipulateArgv.h"

#include "PolyMC.h"
#include "MCMethods/ReplicaExchange.h"
#include "MCMethods/TopologyExchange.h"
#define POLYMC_VERSION 0.720

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

int main(int argc, const char **argcv) {
    std::cout << "PolyMC version: " << POLYMC_VERSION << std::endl;
    arma::arma_version ver;
    std::cout << "ARMA version: "<< ver.as_string() << std::endl;

    std::vector<std::string> argv(argcv, argcv+argc);

    std::string inputfn = "";
    inputfn  = parse_arg(inputfn, "-in",argv);

    // terminate if no input file given
    if (inputfn == "") {
        std::cout << "Error: No input file provided." << std::endl;
        std::exit(0);
    }

    InputRead * input = new InputRead(inputfn);

    if (input->contains_multi("replica_exchange")) {
        ReplicaExchange reex(inputfn,argv);
    }
    else if (input->contains_multi("topology_exchange")) {
        TopologyExchange toex(inputfn,argv);
    }
    else {
        PolyMC polyMC(argv);
        polyMC.simple_execute();
    }

    return 0;
}