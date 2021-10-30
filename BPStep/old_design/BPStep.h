#ifndef __BPSTEP_INCLUDED__
#define __BPSTEP_INCLUDED__

#define PAIR_EVAL_FORWARD 1
/*
Defines whether pair energy evaluations are handled on first or last encounter.
1 for first encounter
*/

#include "../ExtraFuncs.h"
#include "../SO3Methods.h"
#include "IDB.h"

#include "EvalEnergy.h"

#include "EvalEnergy/EE_StiffMat.h"


#include <armadillo> // arma
#include <iostream> // std::cout
#include <string>   // std::string
#include <algorithm>
#include <stdexcept>
#include <vector>

using namespace std;


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

    EvalEnergy* eval_energy_diag;
    std::vector<EvalEnergy*> eval_energy_left;
    std::vector<EvalEnergy*> eval_energy_right;

    arma::mat M0;

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
    double T;
    const double T_ref  = 300;

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
    EvalEnergy * select_EvalEnergy(const std::string & method, const std::vector<double> & params, bool is_diag);

public:
    int            get_id();
    std::string    get_type();
    arma::colvec * get_T0();
    double         get_T0_component(int component);
    arma::mat *    get_R0();

    EvalEnergy *   get_eval_energy_diag();
    EvalEnergy *   get_eval_energy_left(int nl);
    EvalEnergy *   get_eval_energy_right(int nr);

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
    bool        same_EvalEnergy(EvalEnergy * eval1,EvalEnergy * eval2);
//    bool        same_mat(arma::mat & M1, arma::mat & M2);
//    void        set_Ml(int nl, arma::mat & M);
//    void        set_Mr(int nr, arma::mat & M);
//    arma::mat   avg_mat(arma::mat & M1, arma::mat & M2, bool harmonic=false);

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

    void change_T(double new_T);
    void set_T0_subtract(bool subtract=true);

// Accessor
    arma::colvec* get_Theta();
    double        get_Theta_component(int component);
    arma::colvec* get_trial_Theta();
    double        get_energy();
    double        get_trial_energy();

    double get_energy_single(int displ);
    double get_trial_energy_single(int displ);

    double get_torsional_torque();

    bool   local_isotropic_hamiltonian();

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

//protected:
//    double eval_diag_energy(const arma::colvec & Theta);


};

#endif
