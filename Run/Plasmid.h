#ifndef __PLASMID_INDLUCED__
#define __PLASMID_INDLUCED__

// include system classes
#include "../PolyMC.h"
#include "../BPStep/BPStep.h"
#include "../Chain.h"

// include Monte Carlo Steps
#include "../MCStep/MCStep.h"
#include "../MCStep/MCS_Pivot.h"
#include "../MCStep/MCS_CSrot.h"
#include "../MCStep/MCS_CSrot2d.h"
#include "../MCStep/MCS_Slither2d.h"
#include "../MCStep/MCS_CStrans.h"
#include "../MCStep/MCS_Twist.h"

// include Dumps
#include "../dump/Dump.h"
#include "../dump/Dump_Thetas.h"
#include "../dump/Dump_xyz.h"
#include "../dump/Dump_State.h"
#include "../dump/Dump_Stiff.h"
#include "../dump/Dump_PersistenceLength.h"
#include "../dump/Dump_ForceExtension.h"
#include "../dump/Dump_Energy.h"
#include "../dump/Dump_TorsionalStiffness.h"
#include "../dump/Dump_Linkingnumber.h"
#include "../dump/Dump_DistanceMap.h"
#include "../dump/Dump_PlasmidEndpoints.h"
#include "../dump/Dump_HatCurve.h"
#include "../dump/Dump_WritheMap.h"

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

void    run_plasmid(int argc, const char **argv);

#endif
