#ifndef __POLYMC_INCLUDED__
#define __POLYMC_INCLUDED__

//////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////
/*
    Defines to include simulation modes
*/

#define INCLUDE_OPEN_MODE
#define INCLUDE_TWEEZER_MODE
#define INCLUDE_PLECTONEME_MODE
#define INCLUDE_PLASMID_MODE
#define INCLUDE_UMBRELLAPLASMID_MODE
#define INCLUDE_2D_MODES


//////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////


//// include input functions
#include "Input/Argparse.h"
#include "Input/InputRead.h"
#include "Input/InputChoice.h"
#include "Input/GenInputFile.h"

// include seed generator
#include "Extra/Seed.h"

// include system classes
#include "Chain.h"

#ifdef INCLUDE_PLECTONEME_MODE
#include "PlecFreeEnergy/GenPFEConf.h"
#endif


// include Monte Carlo Moves
#include "MCStep/MCStep.h"
#include "MCStep/MCS_ClusterTwist.h"
#include "MCStep/MCS_CSrot.h"
#include "MCStep/MCS_CStrans.h"
#include "MCStep/MCS_PivCon.h"
#include "MCStep/MCS_Pivot.h"
#include "MCStep/MCS_TailTwist.h"
#include "MCStep/MCS_TorquePiv.h"
#include "MCStep/MCS_Twist.h"
#include "MCStep/MCS_zPiv.h"
#include "MCStep/MCS_Kinkxy.h"

// include 2D Monte Carlo Moves
#include "MCStep/MCS_CSrot2d.h"
#include "MCStep/MCS_CSrot_2dpot.h"
#include "MCStep/MCS_Pivot2d.h"
#include "MCStep/MCS_Slither2d.h"

// include Dumps
#include "dump/InitDump.h"

// include Excluded Volume
#include "ExcludedVolume/ExcludedVolume.h"

// include Electrostatics
#include "Electrostatics/Electrostatics.h"
#include "Electrostatics/ESPotential.h"
#include "Electrostatics/Potentials/DebyeHueckel.h"



// include Constraints
#include "Constraints/Constraint.h"
#include "Constraints/Constr_Position.h"
#include "Constraints/Constr_Orientation.h"


// include Unbound
#include "Unbound/Unbound.h"
#include "Unbound/Pair.h"
#include "Unbound//PairStyles/Pair_LJ.h"


// include Geometric and Extra functions
#include "SO3Methods.h"
#include "ExtraFuncs.h"

// include Plectoneme Finder
#include "PlecFinder/PlecFinder.h"

#include <chrono>
#include <armadillo>
#include <cmath>
#include <vector>
//#include <random>

// include electrostatic dump
#include "dump/Dump_ESEnergy.h"


#define POLYMC_MODE_TWEEZER          "tweezer"
#define POLYMC_MODE_TWEEZER          "tweezer"
#define POLYMC_MODE_PLEC             "plec"
#define POLYMC_MODE_PLASMID          "plasmid"
#define POLYMC_MODE_UMBRELLAPLASMID  "umbrellaplasmid"
#define POLYMC_MODE_OPEN             "open"
#define POLYMC_MODE_LIN2D            "lin2d"
#define POLYMC_MODE_CLOSED2D         "closed2d"

#define TWEEZER_SETUP_FREE              "free"
#define TWEEZER_SETUP_LINKFIX           "fix_link"
#define TWEEZER_SETUP_TORSIONAL_TRAP    "torsional_trap"
#define TWEEZER_SETUP_TORQUE            "torque"

#define TWEEZER_SETUP_ID_FREE           0
#define TWEEZER_SETUP_ID_LINKFIX        1
#define TWEEZER_SETUP_ID_TORSIONAL_TRAP 2
#define TWEEZER_SETUP_ID_TORQUE         3

#define TWEEZER_SETUP_LINKFIX_DYNAMICS      "fix_link_dyn"
#define TWEEZER_SETUP_ID_LINKFIX_DYNAMICS   4


#define PLEC_SETUP_FREE                 "free"
#define PLEC_SETUP_LINKFIX              "fix_link"
#define PLEC_SETUP_TORSIONAL_TRAP       "torsional_trap"
#define PLEC_SETUP_TORQUE               "torque"

#define PLEC_SETUP_ID_FREE              0
#define PLEC_SETUP_ID_LINKFIX           1
#define PLEC_SETUP_ID_TORSIONAL_TRAP    2
#define PLEC_SETUP_ID_TORQUE            3

#define LINK_CONSERVATION_VIOLATION_THRESHOLD 1.2

#define DEFAULT_DUMP_DIR "default_output"

//////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////
/*
DEBUG FLAGS
*/
#define POLYMC_TERMINATE_INCONSISTENT_ENERGY
#define POLYMC_TERMINATE_INCONSISTENT_CONFIG

#define SEED_SEQ_LEN 4

//////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////


class PolyMC;

class PolyMC {
//    friend class PolyMC; // not sure what the fuck this is?

protected:
    // Settings
    std::string mode;

/*
    Input Members
*/
    bool        input_given;
    std::string inputfn;
    InputRead * input;

    std::string                 seedstr = "-1";
    std::vector<long long int>  seedseq;
    GenInputFile geninfile;

    std::vector<std::string> argv;
//    int argc;
//    const char **argv;

    std::string restart_file     = "";
    int         restart_snapshot = -1;
    bool        restart_set_link = false;

/*
    Chain Members
*/
    Chain* chain;
    bool chain_initialized = false;

    unsigned num_bp     = 0;
    std::string IDB_fn  = "";
    bool T0_subtract    = false;

    std::string seq_fn  = "";
    std::string seq     = "";

    double temp         = 300;

    double hel_rep_len  = PRESET_HELICAL_REPEAT_LENGTH;
    double sigma        = 0;
    double torque       = 0;
    double force        = 0;
    arma::colvec fdir   = {0,0,1};

    bool subtract_T0    = false;

    bool closed                  = false;
    bool terminus_fixed_tangents = false;
    bool terminus_fixed_triads   = false;

    /*
        In terms of elastic couplings the chain is not closed, but
        pair interactions such as electrostatics consider the chain
        closed
    */
    bool pseudo_closed           = false;

    bool use_cluster_twist = true;
    unsigned num_twist = 0;

/*
    MCStep Members
*/
    std::vector<MCStep*> MCSteps;


/*
    Dump Members
*/
    bool nooutput;
    std::vector<Dump*>   Dumps;
    Dump*                DumpLast;
    std::string dump_dir = DEFAULT_DUMP_DIR;


/*
    Excluded Volume Members
*/
    bool    EV_active           = false;
    bool    EV_repulsion_plane  = false;
    bool    EV_line_closure     = false;

    std::string tweezer_boundary_front = "surface";
    std::string tweezer_boundary_back  = "surface";
    bool bf_surface,bb_surface,bf_line,bb_line;

    ExVol*  EV;
    double  EV_rad = 0;
    bool    EV_check_crossings = true;

/*
    Electrostatics Members
*/
    std::string ElStat_fn   = "";
    bool    ES_active       = false;
    ElStat      * ES;
    ESPotential * ES_Potential;

    double  ES_rho_max          = 0;
    double  ES_rho_min          = 0;
    double  neighbor_skip_dist  = 0;

/*
    Slow Wind
*/
    bool      slow_wind                    = false;
    double    slow_wind_dLK_step           = 0;
    long long slow_wind_steps_per_dLK_step = 10000;
    bool      slow_wind_dump               = false;


/*
    Tweezer Members
*/
    int         tweezer_setup_id = TWEEZER_SETUP_ID_FREE;
    std::string tweezer_setup    = TWEEZER_SETUP_FREE;

/*
    Plec Members
*/
    int         plec_setup_id       = PLEC_SETUP_ID_FREE;
    std::string plec_setup          = PLEC_SETUP_FREE;
    int         plec_num_bp_actual  = 0;
    int         plec_id_first       = 0;
    double      plec_termini_dist   = 0;
    double      plec_loop_frac      = 2./3;

    double      plec_dLK            = 0;
    double      plec_trap_stiffness = 0;

/*
    UmbrellaPlasmid Members
*/
    int         umbrellaplasmid_setup_id       = PLEC_SETUP_ID_FREE;
    std::string umbrellaplasmid_setup          = PLEC_SETUP_FREE;

    double      umbrellaplasmid_dLK            = 0;
    double      umbrellaplasmid_trap_stiffness = 0;

/*
    Lin2d Members
*/
    bool use_slither2d = false;

/*
    General Setup Members
*/

    int  print_every             = 100000;
    bool print_link_info         = false;
    bool print_elstat_info       = false;

    bool check_link              = false;
    int  check_link_every        = 100000;

    arma::cube check_link_backup_triads;
    arma::mat  check_link_backup_bp_pos;
    double     check_link_backup_dLK;

    bool check_consistency             = true;
    long long  check_consistency_every = 1000000;

    bool copy_input = false;

/*
    Run Members
*/
    long long equi  = 0;
    long long steps = 0;

/*
    Print Members
*/
    std::chrono::high_resolution_clock::time_point timer_start;
    std::chrono::high_resolution_clock::time_point timer_finish;
    std::chrono::high_resolution_clock::time_point timer_ref;
    std::chrono::duration<double> timer_elapsed;

/*
    Constraints
*/
    bool constraints_active  = false;
    std::string constraint_fn = "";
    std::vector<Constraint*> constraints;

/*
    Generalized Forces
*/
//    bool GenForce_Active  = false;
    std::string GenForce_fn = "";

/*
    Pair Interactions
*/
    Unbound * unbound;
    bool    pair_interactions_active = false;
    std::string pair_interactions_fn      = "";
    std::string initialize_atom_types_fn  = "";


public:
    PolyMC(const std::vector<std::string> & argv,bool produce_output=true);
    ~PolyMC();

/*
    MC main loop
*/

public:
    bool simple_execute();
    bool run(long long steps,const std::string & run_name,bool dump=true,bool print=true);
    bool run(long long steps,std::vector<MCStep*> run_MCSteps,const std::string & run_name,bool dump=true,bool print=true);

    long long int external_run_steps = 0;
    void init_external_run();
    bool external_run(long long steps,const std::string & run_name,bool dump=true,bool print=true,bool recal_energy=true, bool set_backups=true, bool reset_count=true);
/*
    Getters and Setters
*/
    bool   exchange_configs(PolyMC * other_polymc);
    double get_energy();
    Chain * get_chain();


protected:

    void init_general();
    void init_seed();
    void init_geninfile();
    std::string read_seq(unsigned num_bp);
    void init_dump();
    void copy_input_files();
    void init_ExVol();
    void init_ElStat(const std::string & ElStat_fn);
    bool init_electrostatic_potential(const std::string & type, std::vector<double> & params, double integral_dx, double rho_max);


    void init_constraints(const std::string & constraint_fn, bool append=false);
    void init_constraints();

    void init_GenForce(const std::string & GenForce_fn);

    void init_pair_interactions();

/*
    Preset Protocols
*/

    bool init_open();
    bool init_tweezer();
    bool init_plec();
    bool init_plasmid();
    bool init_lin2d();
    bool lin2d_slither_equi_with_pivot();
    bool init_closed2d();
    bool init_umbrellaplasmid();

/*
    Slow Winds
*/
    bool generic_slow_wind(double dLK_from, double dLK_to, double dLK_step, long long steps_per_dLK_step);
    bool tweezer_slow_wind(double dLK_from, double dLK_to, double dLK_step, long long steps_per_dLK_step);
    bool plec_slow_wind   (double dLK_from, double dLK_to, double dLK_step, long long steps_per_dLK_step);

/*
    Consistency Checks
*/
    bool check_chain_consistency(long long step);
    bool check_link_consistency(long long step);
    void set_link_backup();

/*
    Print Methods
*/
    void print_acceptance_rates(const std::vector<MCStep*> & run_MCSteps);
    std::string time_remaining(int seconds);

    void start_timers();
    void print_state(long long step,long long steps,const std::string & run_name,const std::vector<MCStep*> & run_MCSteps);


//    double run(long long int steps,bool dump=true,bool print=false);

    std::string get_first_word(std::string & str);
    std::vector<std::string> get_str_list(std::string &str);
    void split_string(std::string &str, const std::string &delims, std::vector<std::string> &output);


};

#endif

