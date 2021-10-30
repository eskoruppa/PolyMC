#ifndef __OPEN_ELASTICITY__
#define __OPEN_ELASTICITY__

// include system classes
#include "../PolyMC.h"
#include "../BPStep/BPStep.h"
#include "../Chain.h"

// include Monte Carlo Steps
#include "../MCStep/MCStep.h"
#include "../MCStep/MCS_Pivot.h"
#include "../MCStep/MCS_CSrot.h"
#include "../MCStep/MCS_CStrans.h"
#include "../MCStep/MCS_PivCon.h"
#include "../MCStep/MCS_Twist.h"
#include "../MCStep/MCS_zPiv.h"

// include Dumps
#include "../dump/InitDump.h"
//#include "../dump/Dump.h"
//#include "../dump/Dump_Thetas.h"
//#include "../dump/Dump_xyz.h"
//#include "../dump/Dump_Stiff.h"
//#include "../dump/Dump_PersistenceLength.h"
//#include "../dump/Dump_ForceExtension.h"
//#include "../dump/Dump_Energy.h"
//#include "../dump/Dump_TorsionalStiffness.h"
//#include "../dump/Dump_Linkingnumber.h"
//#include "../dump/Dump_DistanceMap.h"
//#include "../dump/Dump_PlasmidEndpoints.h"
//#include "../dump/Dump_HatCurve.h"

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

void    run_open(int argc, const char **argv);

#endif
