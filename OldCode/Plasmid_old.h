#ifndef __PLASMIDOLD_INCLUDED__
#define __PLASMIDOLD_INCLUDED__

// include system classes
#include "../PolyMC.h"
#include "../BPStep.h"
#include "../Chain.h"
#include "../InteractionDatabase.h"

// include Monte Carlo Steps
#include "../MCStep/MCStep.h"
#include "../MCStep/MCS_Pivot.h"
#include "../MCStep/MCS_CSrot.h"
#include "../MCStep/MCS_CStrans.h"
#include "../MCStep/MCS_Twist.h"
#include "../MCStep/MCS_Reptation_Thermal.h"
#include "../MCStep/MCS_Unwinder.h"
#include "../MCStep/MCS_BranchWinder.h"

// include Dumps
#include "../dump/Dump.h"
#include "../dump/Dump_Thetas.h"
#include "../dump/Dump_xyz.h"
#include "../dump/Dump_Stiff.h"
#include "../dump/Dump_PersistenceLength.h"
#include "../dump/Dump_ForceExtension.h"
#include "../dump/Dump_Energy.h"
#include "../dump/Dump_TorsionalStiffness.h"
#include "../dump/Dump_Linkingnumber.h"
#include "../dump/Dump_DistanceMap.h"
#include "../dump/Dump_PlasmidEndpoints.h"

// include Excluded Volume
#include "../ExcludedVolume/EV_Beads.h"
#include "../ExcludedVolume/ProximityCluster.h"
#include "../ExcludedVolume/ExcludedVolume.h"

#include <chrono>
//using namespace std;
using namespace std::chrono;

#include <armadillo>
#include <cmath>
#include <vector>

#include <omp.h>


void    plasmid();
void    plasmid_replica_exchange();
void    plasmid_unwinding_sampling();

bool    ReplicaExchange(vector<PolyMC*>& Runs, vector<int>& swap_counts);
double  ReplicaExchange_Metropolis(PolyMC* RunA, PolyMC* RunB);
void    ReplicaExchange_Swap(PolyMC* RunA, PolyMC* RunB);

bool RE_run_SERIAL(vector<PolyMC*>& Runs,vector<double>& Temps, long long int steps, long long int sweeps, string& dump_dir);
bool RE_run_OMP(vector<PolyMC*>& Runs,vector<double>& Temps, long long int steps, long long int sweeps, string& dump_dir);


#endif
