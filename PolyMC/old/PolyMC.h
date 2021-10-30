#ifndef __POLYMC_INCLUDED__
#define __POLYMC_INCLUDED__


// include system classes
#include "BPStep/BPStep.h"
#include "Chain.h"

// include Monte Carlo Steps
#include "MCStep/MCStep.h"
#include "MCStep/MCS_Pivot.h"
#include "MCStep/MCS_CSrot.h"
#include "MCStep/MCS_CStrans.h"
#include "MCStep/MCS_Twist.h"
#include "MCStep/MCS_Unwinder.h"

// include Dumps
#include "dump/Dump.h"
#include "dump/Dump_Thetas.h"
#include "dump/Dump_xyz.h"
#include "dump/Dump_Stiff.h"
#include "dump/Dump_PersistenceLength.h"
#include "dump/Dump_ForceExtension.h"
#include "dump/Dump_Energy.h"
#include "dump/Dump_TorsionalStiffness.h"
#include "dump/Dump_Linkingnumber.h"
#include "dump/Dump_DistanceMap.h"
#include "dump/Dump_PlasmidEndpoints.h"

// include Excluded Volume
#include "ExcludedVolume/EV_Beads.h"
#include "ExcludedVolume/ProximityCluster.h"
#include "ExcludedVolume/ExcludedVolume.h"

#include <chrono>
//using namespace std;
using namespace std::chrono;

#include <armadillo>
#include <cmath>
#include <vector>

class PolyMC;

class PolyMC {

    protected:
    long long int step;


    bool setup_finalized=false;
    bool kill_if_EV_voilated=true;
    int  cal_Lk_every = 10000;

    // Chain
    Chain* chain;

    // MCSteps
    vector<MCStep*> MCSteps;

    //Dumps
    vector<Dump*>   Dumps;

    // Excluded Volume
    protected:
    bool    EV_active=false;
    ExVol*  EV;


    // Settings
    string type;


    // Variables




public:
    PolyMC();
    ~PolyMC();

    void setup(string type, string interaction_file, int num_bp, string sequence, double sigma, double T);
    void set_ExVol(double EV_rad, bool repulsion_plane);
    void set_dumps(string dir, int dump_every);

/*
    Mutators
*/
public:
    void set_T(double T);

protected:
    void finalize_setup();

public:

    Chain* get_chain();
    ExVol* get_EV();


    double run(long long int steps,bool dump=true,bool print=false);


protected:
    double  sigma_unwound;
    int     steps_unwinding;
    int     steps_winding;
    int     step_equi;
    int     n_split_winding;

    MCS_Unwinder * MCSU;

public:
    void   unwinding_setup(ExVol * EV, double sigma_unwound, int steps_unwinding, int steps_winding, int step_equi, int n_split_winding);
    double unwinding_run();




};

#endif
