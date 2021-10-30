#ifndef __SIM2D_INDLUCED__
#define __SIM2D_INDLUCED__

// include system classes
#include "../PolyMC.h"
#include "../BPStep/BPStep.h"
#include "../Chain.h"

// include Monte Carlo Steps
#include "../MCStep/MCStep.h"
#include "../MCStep/MCS_Pivot.h"
#include "../MCStep/MCS_CSrot.h"
#include "../MCStep/MCS_CSrot_2dpot.h"
#include "../MCStep/MCS_CSrot2d.h"
#include "../MCStep/MCS_Pivot2d.h"
#include "../MCStep/MCS_Slither2d.h"
#include "../MCStep/MCS_CStrans.h"
#include "../MCStep/MCS_Twist.h"

// include Dumps
#include "../dump/InitDump.h"

// include Excluded Volume
#include "../ExcludedVolume/EV_Beads.h"
#include "../ExcludedVolume/ProximityCluster.h"
#include "../ExcludedVolume/ExcludedVolume.h"

// include input functions
#include "../Input/Argparse.h"

#include <chrono>
//using namespace std;
using namespace std::chrono;

#include <armadillo>
#include <cmath>
#include <vector>

#include <omp.h>

void    run_teardrop2d(int argc, const char **argv);
void    run_linear2d(int argc, const char **argv);

#endif
