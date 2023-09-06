#ifndef __BPSTEP_INCLUDED__
#define __BPSTEP_INCLUDED__

#include "../ExtraFuncs.h"
#include "../SO3Methods.h"
#include "IDB.h"

#include "EvalEnergy.h"
#include "EvalEnergy/EE_StiffMat.h"
#include "EvalEnergy/EE_KinkXY.h"


#include <armadillo> // arma
#include <iostream> // std::cout
#include <string>   // std::string
#include <algorithm>
#include <stdexcept>
#include <vector>

///////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////
#define PAIR_EVAL_FORWARD 1
/*
Defines whether pair energy evaluations are handled on first or last encounter.
1 for first encounter
*/



#define REF_TEMP 300
#define REF_KBT  4.114

#define STIFFMAT_FACTOR          0.5

#define BPS_MC_COV_DEFAULT_STEPS 1e6
#define BPS_MC_COV_DEFAULT_SIGMA 0.1

#define BPS_ISOTROPIC_BENDING_NUM_SAMPLED 1e5
#define BPS_ISOTROPIC_BENDING_EQUAL_THRES 1e-8

///////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////


#define BPS_USE_EVALENERGY
#define INCLUDE_GENERALIZED_FORCE

///////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////



//using namespace std;


class BPStep;

class BPStep
{
// member variables
protected:
    int          id;     // id in chain
    std::string  type;     // basepair type e.g. AT ; refers to sequence on reference strand
    bool         closed;
    bool         local_isotropic;

    std::vector<BPStep*>  BPSl;  // Pointers to left neighbour BPSteps
    std::vector<BPStep*>  BPSr;  // Pointers to right neighbour BPSteps

    bool T0_subtract; // by default the intrinsic rotational components are accounted for by a left hand rotation matrix (R0)

    arma::colvec* T0;     // Theta0 - intrinsic rotations
    arma::mat*    R0;
    arma::mat*    R0_T;

    arma::colvec  T0_store;     // Theta0 - intrinsic rotations
    arma::mat     R0_store;
    arma::mat     R0_T_store;

    arma::mat     cov;

    double eval(const arma::colvec & Theta);


//    arma::mat*    M0;     // M_onsite - local stiffness matrix
//    arma::mat     M0_store;     // M_onsite - local stiffness matrix
//
//    std::vector<arma::mat*>  Ml;         // M_lower  - pointers to stiffness matrices coupling to BPS with lower index
//    std::vector<arma::mat*>  Mr;         // M_upper  - pointers to stiffness matrices coupling to BPS with larger index
//
//    std::vector<arma::mat>   Ml_store;    // M_lower  - stiffness matrices coupling to BPS with lower index
//    std::vector<arma::mat>   Mr_store;    // M_upper  - stiffness matrices coupling to BPS with larger index

    int coupling_range;  // range of coupling matrices (equals amount matrices contained in Ml and Mu)
    int coup_left;
    int coup_right;

    double disc_len;
    double temp;
    const double ref_temp  = 300;

    double kT = 4.114;
    const double ref_kT = 4.114;

    double               energy_0;
    std::vector<double*> energy_l;  // Pointers to current energetic contribution stemming from the interaction with certain left neighbor
    std::vector<double*> energy_r;  // Pointers to current energetic contribution stemming from the interaction with certain right neighbor

    double               trial_energy_0;
    std::vector<double*> trial_energy_l;
    std::vector<double*> trial_energy_r;

    arma::colvec    Theta;
    arma::colvec    trial_Theta;
    bool            trial_pending;
    bool            trial_eval_pending;

    #ifdef INCLUDE_GENERALIZED_FORCE
    bool            generalized_force_active;
    arma::colvec    generalized_force;

    double          generalied_force_energy;
    double          generalied_force_trial_energy;
    #endif


// member functions
public:
    BPStep(int id_in_chain, const std::string& full_seq, IDB& database, double disc_len, double temp, bool closed);
    ~BPStep();

public:
    bool init_neighbors(std::vector<BPStep*>& BPS_all);
    void check_neighbors();
protected:
    std::string find_type(int id_in_chain, const std::string& full_seq, int coup_range);
    bool   set_params(const std::string& type, int coup_range, IDB& data);
    arma::mat vec2mat(std::vector<double> params);


public:
    int            get_id();
    std::string    get_type();
    arma::colvec * get_T0();
    double         get_T0_component(int component);
    arma::mat *    get_R0();

    arma::mat get_cov();
//    arma::mat *    get_M0();
//    arma::mat *    get_Mr(int nr);
//    arma::mat *    get_Ml(int nl);
//    arma::mat *    get_Mc(intdispl);

    int get_range();
    int get_coup_left();
    int get_coup_right();

/////////////////////////////////
// Coupling Matrix consistency //

public:
    bool check_coupling_consistency(bool avg_conflicts, bool harmonic=false);

    void print_M_ptr();

protected:
//    bool        same_mat(arma::mat & M1, arma::mat & M2);
//    void        set_Ml(int nl, arma::mat & M);
//    void        set_Mr(int nr, arma::mat & M);
//    arma::mat   avg_mat(arma::mat & M1, arma::mat & M2, bool harmonic=false);


    bool                equal_params(const std::vector<double> & params1,const std::vector<double> & params2);
    std::vector<double> avg_params(const std::vector<double> & params1,const std::vector<double> & params2, bool harmonic=false);

////////////////////////////////////////////////////////////////////////////////////
// Energy Eval Functions ///////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////


// Main Eval Methods

/*
    The protocol for energy evaluation is to call the three trial functions
    for all BPSteps whose internal defrees of freedom changed upon the trial
    move. Those functions have to be called in consecutive order and
    eval_delta_energy cannot be called for any BPStep before the propose_move
    was executed for all involved BPStep. The same hold for set_move and
    eval_delta_energy.

    propose_move (trial_Theta)
    -   Sets the new degrees of freedom as a trial parameters. The new value
        is only stored as a trial and can only become the actual value
        if the move is accepted. In the BPSteps for which this function was
        called a flag will be set indicating to the other BPSteps that
        that this instance will evaluate an energy during the current trial
        round. This is important when coupling terms are involved to prevent
        double evaluating of those coupling terms

    eval_delta_energy ()
    -   Evaluates the energy of the trial move and stores the trial energies.
        Returns the change of energy with respect to the previous configuration.
        However, the return value of this function will not always be the complete
        change of energy associated with the repective instance. In order to
        avoid double counting of coupling energies, those terms will only be considered
        if the coupling partner either already did evaluate its energy at the current
        trial run or was not flagged to be evaluated. Therefore this function will
        not be reliable for individual instances but the commulative energy difference
        will be correct.
        If the total energetic contribution of a single instance needs to be accessed the
        function BPStep::get_trial_energy_single() can be utilized.

    set_move (true/false)
        Has to be called both in the case of acceptance and rejectance to finalize the
        pending trial move. Resets the flags and confirms the trial Theta and energies to
        become the current values.
*/
public:

// Core

//    void propose_move(arma::colvec& Theta);
    void   propose_move(const arma::mat& T1, const arma::mat& T2);
    double eval_delta_energy();
    double eval_new_energy();
    void   set_move(bool accepted);

    bool eval_pending();
    bool pending();

// Extra

    bool pending_closed();
    void close_pending();

// Mutator

    void change_T(double new_temp);
    void set_T0_subtract(bool subtract=true);

    #ifdef INCLUDE_GENERALIZED_FORCE
    void   set_generalized_force(const arma::colvec & gen_force);
    void   deactivate_generalized_force();
    double eval_generalized_force(const arma::colvec & Theta);
    double eval_change_generalized_force(const arma::colvec & Theta);
    #endif


    bool twist_energy_active = true;
    void deactivate_twist_energy();
    void reactivate_twist_energy();

// Accessor
    arma::colvec* get_Theta();
    arma::colvec  get_FullTheta();
    double        get_Theta_component(int component);
    arma::colvec* get_trial_Theta();
    double        get_energy();
    double        get_trial_energy();

    double get_energy_single(int displ);
    double get_trial_energy_single(int displ);

    double get_torsional_torque();

    bool   local_isotropic_hamiltonian();

    std::vector<double> get_status_diag();
    std::string *       get_method_diag();

    EvalEnergy * get_EvalEnergy_diag();

    // methods to relate deformation vector Theta (excess rotation) to triads for given BPS
    arma::colvec Triads2Theta       (const arma::mat & T1,const arma::mat & T2);
    arma::colvec Triads2FullTheta   (const arma::mat & T1,const arma::mat & T2);
    arma::mat    Theta2Rotmat       (const arma::colvec & Theta);
    arma::mat    FullTheta2Rotmat   (const arma::colvec & FullTheta);

//////////////////////////////////////
////// Extract Current Energies //////

public:
    void    energy_extract_select();
    bool    energy_extract_is_selected();
    double  energy_extract();

protected:
    bool    energy_extract_selected;


//////////////////////////////////////
//////////////////////////////////////


// Init

/*
    These methods should only be used during the initialization stage when the chain
    is build. They facilitate the mutual sharing of coupling energies.
*/
public:
    void wire_energies();
protected:
    void assign_energy_l_pointers(int id, double* el, double* pel);
    void assign_energy_r_pointers(int id, double* er, double* per);



#ifdef BPS_USE_EVALENERGY
///////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////
////////////////////// EVAL ENERGY ////////////////////////////////////////////////////////


protected:
    EvalEnergy* eval_energy_diag;
    std::vector<EvalEnergy*> eval_energy_left;
    std::vector<EvalEnergy*> eval_energy_right;

    EvalEnergy * select_EvalEnergy(const std::string & method, const std::vector<double> & params, bool is_diag);
    bool         same_EvalEnergy(EvalEnergy * eval1,EvalEnergy * eval2);

public:
    EvalEnergy *   get_eval_energy_diag();
    EvalEnergy *   get_eval_energy_left(int nl);
    EvalEnergy *   get_eval_energy_right(int nr);

#endif

    std::vector<double> get_eval_energy_params(int id);


///////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////
/////////////////// INTERNAL ENERGY EVAL //////////////////////////////////////////////////
#ifndef BPS_USE_EVALENERGY

public:

    std::vector<double> params_diag;
    arma::mat           params_diag_mat;
    unsigned            coup_diag_method;

    std::vector<std::vector<double>>    params_left;
    arma::cube                          params_left_mat;
    std::vector<unsigned>               coup_left_method;

    std::vector<std::vector<double>>    params_right;
    arma::cube                          params_right_mat;
    std::vector<unsigned>               coup_right_method;


protected:
    unsigned            get_method_id(const std::string & method);
    std::vector<double> set_params(unsigned method, std::vector<double> params);
    bool                allow_params_averaging(unsigned method);
    std::vector<double> set_params_temp(double new_temp,unsigned method, std::vector<double> params);
    bool                check_isotropic_bending();
    bool                cal_MC_isotropic_bending();
    arma::mat           cal_cov();
    arma::mat           cal_MC_cov(long long int steps=BPS_MC_COV_DEFAULT_STEPS, double sigma=BPS_MC_COV_DEFAULT_SIGMA);

    double              cal_beta_energy_diag(const arma::colvec & Theta);
    double              cal_beta_energy_left(const arma::colvec & Theta1,const arma::colvec & Theta2, unsigned id_left);
    double              cal_beta_energy_right(const arma::colvec & Theta1,const arma::colvec & Theta2, unsigned id_right);






    /*
        Stiffmat Methods
    */

    std::vector<double> stiffmat_set_params(std::vector<double> params);
    std::vector<double> stiffmat_set_params_temp(double new_temp,std::vector<double> params);

    double stiffmat_cal_beta_energy_diag(const arma::colvec & Theta);
    double stiffmat_cal_beta_energy_offdiag(const arma::colvec & Theta1,const arma::colvec & Theta2, arma::mat M);

#endif


};

#endif
