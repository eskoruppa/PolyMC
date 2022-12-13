#ifndef __STATE_INCLUDED__
#define __STATE_INCLUDED__

#include "SO3Methods.h"
#include "ExtraFuncs.h"
#include "BPStep/BPStep.h"
#include "BPStep/IDB.h"

////////////////////////////////////////////
// include Load Restart
#include "Input/WriteLoadRestart.h"


#define _USE_MATH_DEFINES
#define ONE_OVER_TWO_PI     0.15915494309189533576888376337251436

#include <armadillo> // armadillo
#include <iostream>  // input output
#include <string>    // string
#include <algorithm>
#include <stdexcept> // exceptions
#include <fstream>   // read file
#include <cmath>     // math lib
#include <math.h>
#include <vector>    // vector


#define MAX_POSITION_DISCREPANCY            1e-8
#define MAX_SO3_DETERMINANT_DISCREPANCY     1e-10
#define AVG_COUP_MISSMATCH          true
#define AVG_COUP_HARMONIC           false

#define PRESET_HELICAL_REPEAT_LENGTH 3.57  // 10.5 basepairs (0.34*10.5)

#define CHAIN_DEBUG


class Chain;

class Chain {

//////////////////////
// member variables //
protected:

    // Energy related members
    const double T_ref  = 300;
    const double kT_ref = 4.114;
    double kT           = 4.114;
    double beta         = 1./4.114;
    double invT_fac     = 1;
    double T            = 300;

    // Preset average intrinsic twist density
    double preset_intrinsic_twist_density = 36./180.*M_PI/0.34;
    double avg_intrinsic_twist_density    = preset_intrinsic_twist_density;

    // Preset Helical Repeat
    double preset_helical_repeat_length = PRESET_HELICAL_REPEAT_LENGTH;
    double helical_repeat_length        = preset_helical_repeat_length;

    // Interaction details
    IDB data;
    std::vector<BPStep*>  BPS;
    unsigned interaction_range;
    bool     T0_subtract = false;
    // 0 for pure onsite couplings.
    // specified by Interaction file, unless
    // overwritten explicitely with
    // set_interaction_range(unsigned int ir)


///////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////
/*
    members related to linking number and torque
*/

protected:
    // Current Linking Number
    double  dLK = 0;

    double  dLK_propose = 0;
    arma::mat terminus_ref_triad;

public:

    void        allow_free_endbead_rotation();
    double      get_dLK();
    arma::mat*  get_terminus_ref_triad();

    void        set_dLK(double dLK);
    void        set_terminus_ref_triad(const arma::mat & ref_triad);

protected:
    // Options
    bool link_torsional_trap = false;
    bool link_const_torque   = false;
    bool link_constrained    = false;
    bool link_track_terminus = false;

public:
    bool torsional_trap_active();
    bool const_torque_active();
    bool track_terminus_active();

    void fix_link(bool fix=true);
    bool link_fixed();

protected:
    // Torsional Trap
    double    torstrap_dLK_fix = 0;
    double    torstrap_stiff   = 0;
    double    torstrap_dLK_aim = 0;

    double    torstrap_torque_sum;
    long long torstrap_torque_count;

///////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////

public:
    void   impose_torsional_trap(double dLK_fix, double trapstiff=-1);
    double get_torstrap_dLK_fix();
    void   set_torstrap_dLK_fix(double dLK);
    double get_torstrap_dLK_aim();
    double get_torstrap_trapstiff();

    double measure_current_torque();
    double get_torque_measure_accu(bool reset_sum=true);


protected:
    // Constant Torque
    double  torque              = 0;
    double  beta_torque         = 0;

public:
    void impose_torque(double torque);
    void change_torque(double torque);


////// Terminus changing Move ///////
protected:
//    bool   termrot_pending = false;
    double termrot_rotation_dLK_new;
    double termrot_torque;
    double termrot_torque_new;

public:
    double propose_terminus_rotation(const arma::mat & T_newlast, double change_angle);
    void   set_terminus_rotation(bool accept);

    void termrot_set_T(double new_T);



///////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////


////    double  trap_stiff         = 0;
//
//    bool    link_constrained   = false;
//
//
//    bool    torque_fixed       = false;
//    bool    torque_measure     = false;
//
//
////    double  dLK_current        = 0;
//    double  torque_trap_dLKaim = 0;
//    arma::mat ref_last_triad;
//
//public:
//
//    // Fixed Linking Number
//    void        fix_link(double dLK);
//    bool        link_fixed();
//
//
//
//    double      get_torque_trap_stiff();
////    double      get_dLK();
//    double      set_dLK(double dLK);
////    double      get_dLK_fix();
//    void        set_dLK_fix(double dLK);
//    arma::mat*  get_ref_last_triad();
//
//
//
//    bool        torque_fix_active();
//    bool        torque_measure_active();
//    double      get_fixed_torque();
//    double      measure_torque();
//    double      get_torque_trap_dLKaim();
//
//    void set_fix_torque(double tor);
//    void set_measure_torque(double C=100, double trapstiff=0);
//    void remove_torque();
////    void change_torque(double tor);
//
//protected:
//    // Torque Statistics
//    double        torque_measure_sum;
//    long long int torque_measure_count;
//    double        torque_measure_partial_sum;
//    long long int torque_measure_partial_count;
//
//public:
//    double torque_measure_add(double tor);
//    double torque_measure_read();
//    double torque_measure_read_partial(bool reset_sum=true);
//    double get_torsional_torque();



///////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////


protected:
    // Force Constraint
    double  force_constrained   = false;
    double  force               = 0;
    arma::colvec force_dir      = {0,0,1};
    arma::colvec beta_force_vec = {0,0,1};

    // Termini Constraint
    bool termini_fixed           = false;
    bool termini_fix_radial      = false;
    bool last_fixed_orientation  = false;
    bool first_fixed_orientation = false;

    // chain length
    unsigned num_bp;
    unsigned num_bps;

    // topology
    bool  closed_topology        = false;

    /*
        pseudo_closed_topology:
            In terms of elastic couplings the chain is not closed, but
            pair interactions such as electrostatics consider the chain
            closed
    */
    bool  pseudo_closed_topology = false;

    std::string config_type = "unspecified";
    /*
     * closed topology controls whether the last element is connected to the
     * first element. If true energy evaluation as well as MC moves are allowed
     * across the termini
     * Always set this variable with the function set_topology(bool closed)
     * to automatically account for  all related settings
     */

    // principle axis; only relevant if tolology is not closed
    arma::colvec  principle_axis = {0,0,1};

    // discertization length
    double  disc_len = 0;
    double  contour_len = 0;

    bool conf_initialized;
    // monomer sequence
    std::string  seq;

    // Main production run configuration
    // MC steps are performed on this configuration
    arma::mat  bp_pos;
    arma::cube triads;
    //arma::mat  Omegas;
    //arma::mat  energy;


    /*
    TODO: These backup states could be deleted. They should rather be stored in
    the dedicated Excluded Volume Objects.
    */
    // Backup configuration for excluded volume checks
    // If a excluded volume violation is detected the main
    // configuration is replaced with this checkpoint
    // Will only be initialized if EVcheck is active.
    arma::mat  bp_pos_EVcheck;
    arma::cube triads_EVcheck;
    //arma::mat  Omegas_EVcheck;
    //arma::mat  energy_EVcheck;
    bool EVcheck;




///////////////////////// LINK CHECK /////////////////////////////////
    // Backup configuration for linking number consistency check
    // If a linking number violation is detected the main
    // configuration is replaced with this checkpoint
    // Will only be initialized if linkcheck is active.
protected:
    bool linkcheck;
    arma::mat  bp_pos_linkcheck;
    arma::cube triads_linkcheck;
    //arma::mat  Omegas_linkcheck;
    //arma::mat  energy_linkcheck;


protected:
    /*
        Average stiffness matrix and average covariance matrix
    */
    arma::mat avg_stiff;
    arma::mat avg_cov;
    arma::mat avg_chol;


//////////////////////
// member functions //

// ACCESSOR METHODS
public:
    std::vector<BPStep*>* get_BPS();
    arma::mat*  get_bp_pos();
    arma::cube* get_triads();

    double get_disc_len();
    double get_contour_len();

    int get_num_bp();
    int get_num_bps();

    bool local_isotropic_hamiltonian();

    double      get_kT();
    double      get_beta();
    double      get_T();
    double      get_T_ref();
    unsigned    get_interaction_range();

    void        cal_avg_stiff(bool inittrue=false);
    arma::mat   get_avg_stiff();
    arma::mat   get_avg_cov();
    arma::mat   get_avg_chol();

    bool         force_active();
    double       get_force();
    arma::colvec get_force_dir();
    arma::colvec get_beta_force_vec();

    bool        fixed_termini();
    bool        fixed_termini_radial();
    bool        fixed_first_orientation();
    bool        fixed_last_orientation();
    bool        topology_closed();
    bool        topology_pseudo_closed();

    std::string get_config_type();
    std::string get_sequence();

//////////////////////////////////////////////////////////
//////////// CALCULATE AND ACCESS ENERGIES ///////////////

public:
    void    cal_energy_propose(int from, int to);
    double  cal_energy_eval   (int from, int to);
    double  cal_energy         (int from, int to); //this move is a combination of the propose and eval move
    void    cal_energy_set_move(int from, int to,bool accept);

    void    set_new_energy();

    double  extract_energy(int from, int to);
    double  extract_energy();

    double  extract_true_energy();

protected:
    void    extract_energy_select(int from, int to);

public:
    bool    check_energy_consistency();
    double  recal_energy();


// CHAIN CONSISTENCY METHODS

public:
    bool config_consistent();
    void restore_consistency();

// MUTATOR METHODS
public:
    Chain(std::string interaction_file);
    ~Chain();

    void set_helical_repeat_length(double hel_rep);
    void set_intrinsic_twist_density(double twist_density);

    void set_T(double temp);
    void set_T0_subtract(bool subtract);
    void set_interaction_range(unsigned ir);

    void set_force (double f, const arma::colvec& dir={0,0,1});

//    void fix_link  (bool fix=true);

    void fix_termini_position(bool fix_radial=false, bool fix=true);
    void fix_termini_orientation(bool fix=true);

protected:
    void set_closed_topology(bool closed);

public:
    void set_pseudo_closed_topology(bool closed = true);


// CONF GENERATORS
public:
    bool gen_linear  (unsigned number_bp, const std::string& sequence,double supercoiling_density=0, const arma::colvec direction={0,0,1});
    bool gen_circular(unsigned number_bp, double supercoiling_density, const std::string& sequence, bool closed_topology = true);
    bool gen_conf_from_restart(std::string restart_fn, int snapshot=-1 ,std::string required_config_type="");
    bool init_custom_open_conf(arma::mat pos, arma::cube triads, const std::string& sequence);

    bool set_orientation_first_triad(arma::colvec& tangent);
    bool set_orientation_last_triad(arma::colvec& tangent);

    bool set_config(arma::mat* bp_pos, arma::cube* triads,bool closed,bool init=false);

    void set_Delta_Lk(double new_dLK);
    void set_Delta_Lk(double new_dLK, int id_first, int id_last=-1);
    void set_sigma(double sigma);
    void set_sigma(double sigma, int id_first, int id_last=-1);


    double sigma2dLk(double sigma);
protected:
    void init_BPS();
    void init_energies();
    void set_seq(const std::string& sequence, int num_bp);



////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////  OBSERVABLES ///////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////

public:

    bool check_link_conservation();

    double cal_quick_writhe();
    double cal_fuller_writhe(int from, int to, arma::colvec& z_dir);
    double cal_twist(int from, int to);
    double cal_twist();

    double gauss_writhe(double density);
    double cal_langowski_writhe_1a();
    double cal_langowski_writhe_1a(int front_extend, int back_extend);
    double cal_langowski_writhe_1a(double density);
    double cal_langowski_writhe_1a(const arma::mat& pos,bool closed);

    void langowski_writhe_elements(arma::mat* writhe_elem);
    void langowski_writhe_elements(arma::mat* writhe_elem, bool closed);
    void langowski_writhe_elements(arma::mat* writhe_elem, double seg_size);

    double rel_extension(arma::colvec& z_dir);
    double extension(arma::colvec& z_dir);

    arma::mat distmap(double density);

    arma::colvec cal_static_persistence_length(int from, int to, int m_max);


};

#endif

