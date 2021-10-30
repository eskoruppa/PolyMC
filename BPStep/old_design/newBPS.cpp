#include "BPStep.h"
/*
 *  This Class stores the information about a certain base pair step.
 *
 *  Contains informations about groundstate deformation Theta0. Access by pointer to arma::colvec
 *      access T0:
 *          - get_T0()
 *
 *  Other members:
 *      int id                      - id along chain
 *      string type                 - base pair step flavor stored in lower case (eg 'at')
 *      int coupling_range          - range of off-site couplings
 *      bool coupling_consistent    - checks whether forward and backward couplings (upper and lower)
 *                                    are consistent
 *
 */


BPStep::BPStep(int id_in_chain, const std::string& full_seq, IDB& database, double disc_len, double temp ,bool closed)
: id(id_in_chain), closed(closed), trial_pending(false), disc_len(disc_len), temp(temp), T0_subtract(false)
{
    coupling_range = database.get_interaction_range();
    type = find_type(id, full_seq, coupling_range);
    kBT      = REF_KBT;
    ref_temp = REF_TEMP;
    set_params(type, coupling_range, database);
    Theta = arma::zeros(3);
    
    change_T(temp)
    cov   = cal_cov();

}

BPStep::~BPStep(){

}

////////////////////////////////////////////////////////////////////////////
/////////// ACCESSOR METHODS/ //////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////

int  BPStep::get_id() {
    return id;
}

std::string BPStep::get_type() {
    return type;
}

arma::colvec * BPStep::get_T0() {
    return T0;
}

double BPStep::get_T0_component(int component) {
    return (*T0)(component);
}

arma::mat * BPStep::get_R0() {
    return R0;
}

arma::mat BPStep::get_cov() {
    #if BPS_USE_EVALENERGY == 1
    return eval_energy_diag->get_cov();
    #else
    return cov;
    #endif
}

#if BPS_USE_EVALENERGY == 1
EvalEnergy * BPStep::get_eval_energy_diag() {
    return eval_energy_diag;
}

EvalEnergy * BPStep::get_eval_energy_left(int nl) {
    if (nl < coup_left) {
        return eval_energy_left[nl];
    }
    throw std::logic_error("Requesting BPStep::get_eval_energy_left left coupling index out of range!");
}

EvalEnergy * BPStep::get_eval_energy_right(int nr) {
    if (nr < coup_right) {
        return eval_energy_right[nr];
    }
    throw std::logic_error("Requesting BPStep::get_eval_energy_right right coupling index out of range!");
}
#endif

int  BPStep::get_range() {
    return coupling_range;
}

int  BPStep::get_coup_left() {
    return coup_left;
}

int  BPStep::get_coup_right() {
    return coup_right;
}


////////////////////////////////////////////////////////////////////////////
/////////// MUTATOR METHODS/ ///////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////


void BPStep::change_T(double new_temp) {

    #if BPS_USE_EVALENERGY == 1
    eval_energy_diag->set_temp(Tnew);
    for (int i=0;i<coup_left;++i) {
        eval_energy_left[i]->set_temp(Tnew);
    }
    for (int i=0;i<coup_right;++i) {
        eval_energy_right[i]->set_temp(Tnew);
    }
    #else
    params_diag_mat = params_diag_mat*temp/new_temp;
    for (unsigned i=0;i<coup_left;i++) {
        params_left_mat.slice(i) *= new_temp/temp;
    }
    for (unsigned i=0;i<coup_right;i++) {
        params_right_mat.slice(i) *= new_temp/temp;
    }

    params_diag = set_param_temp(new_temp,coup_diag_method, params_diag);
    for (unsigned i=0;i<coup_left;i++) {
        params_left[i] = set_param_temp(new_temp,coup_left_method[i], params_left[i]);
    }
    for (unsigned i=0;i<coup_right;i++) {
        params_right[i] = set_param_temp(new_temp,coup_right_method[i], params_right[i]);
    }
    #endif

    temp = new_temp;
    kBT  = REF_KBT*temp/REF_TEMP;

    #if BPS_USE_EVALENERGY == 1
    cov = eval_energy_diag->get_cov();
    #else
    cov   = cal_cov();
    #endif
}

void BPStep::set_T0_subtract(bool subtract) {
    T0_subtract = subtract;
    if (T0_store(0)==0 && T0_store(1)==0 && T0_store(2)==0) {
        T0_subtract=true;
    }
}


////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////
///////////// INIT METHODS /////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////

std::string BPStep::find_type(int bps, const std::string& full_seq, int coup_range) {
    int         num_bp = int(full_seq.length());
    int         num_bps;
    int         shift;
    int         seglen = (1+coup_range)*2;
    std::string assign_type;

    bool print_map=false;

    // Consistency Check and Output //
    if (num_bp != (int)full_seq.length()) {
        throw std::domain_error("BPStep::find_type(): Sequence length does not coincide with number of basepairs.");
    }
    if (seglen > num_bp) {
        throw std::domain_error("BPStep::find_type(): Interaction length too long for given full_seq.");
    }

    ////////////////////////////////////////////////////////////////////////////
    ////////// SETTING BPStep //////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////////

    if (!closed) {
    /*
    --- OPEN ---
    If the topology is open variable segemts x are added at the endpoints to use averaged
    interactions.
    */
        num_bps=int(num_bp-1);
        ////////////////////////////////
        //////// DETERMINE TYPE ////////
        ////////////////////////////////
        // termini
        if ( int(bps-coup_range) < 0 || int(bps+coup_range) >= num_bps ) {
            // close to left terminus
            if (int(bps-coup_range) < 0) {
                shift = coup_range-bps;
                assign_type = "";
                for (int i=0;i<shift;++i) {
                    assign_type=assign_type+"x";
                }
                assign_type = assign_type+full_seq.substr(0,seglen-shift);
                if (print_map) {std::cout << bps << " " << assign_type << "   <--" << shift << std::endl;}
            }
            // close to right terminus
            if (int(bps+coup_range) >= num_bps) {
                shift = bps+coup_range+1-num_bps;
                assign_type = full_seq.substr(bps-coup_range,seglen-shift);
                for (int i=0;i<shift;++i) {
                    assign_type=assign_type+"x";
                }
                if (print_map) {std::cout << bps << " " << assign_type << "   <--" << shift << std::endl;}
            }
        }
        // central
        else {
            assign_type = full_seq.substr(bps-coup_range,seglen);
            if (print_map) {std::cout << bps << " " << assign_type << std::endl;}
        }

    }
    else {
    /*
    --- CLOSED ---
    If the topology is closed the endpoints will be periodically connected.
    */
        num_bps=(int)(num_bp);

        ////////////////////////////////
        //////// DETERMINE TYPE ////////
        ////////////////////////////////
        // termini
        if ( int(bps-coup_range) < 0 || int(bps+coup_range) >= (num_bps-1) ) {
            // close to left terminus
            if (int(bps-coup_range) < 0) {
                shift = coup_range-bps;
                assign_type = full_seq.substr(num_bps-shift,shift);
                assign_type = assign_type+full_seq.substr(0,seglen-shift);
                if (print_map) {std::cout << bps << " " << assign_type << "   <--" << shift << std::endl;}
            }
            // close to right terminus
            if (int(bps+coup_range) >= (num_bps-1)) {
                shift = bps+coup_range+2-num_bps;
                assign_type = full_seq.substr(bps-coup_range,seglen-shift);
                assign_type = assign_type+full_seq.substr(0,shift);
                if (print_map) {std::cout << bps << " " << assign_type << "   <--" << shift << std::endl;}
            }
        }
        // central
        else {
            assign_type = full_seq.substr(bps-coup_range,seglen);
            if (print_map) {std::cout << bps << " " << assign_type << std::endl;}
        }
    }
    return assign_type;
}

////////////////////////////////////////////////////////////////////////////
///////////////// ASSIGN INTERACTION PARAMETERS ////////////////////////////

bool BPStep::set_params(const std::string& type, int coup_range, IDB& data) {

    int coup = int(type.length()/2-1);

    if (coup != coup_range) {
        throw std::domain_error("BPStep::set_params(): coupling range does not agree with type sequence length");
        return false;
    }

    coup_left  = coup_range;
    coup_right = coup_range;

    // Set left and right coupling ranges and ascertain whether the step is a terminus step
    bool terminus = false;
    for (int i=0;i<coup_range;i++) {
        if (type[i]=='x') {
            --coup_left;
            terminus = true;
        }
        if (type[type.length()-1-i]=='x') {
            --coup_right;
            terminus = true;
        }
    }

    // Retrieve IDBUnit from data and assign member stiffness members
    IDBU * idbu = data.get_IDBU(type);

    /*
        Init energy evaluation classes derived from EvalEnergy
    */

    #if BPS_USE_EVALENERGY == 1
    eval_energy_diag = select_EvalEnergy(idbu->methods[coup_range],idbu->interactions[coup_range],true);
    for (int i=0;i<coup_left;++i) {
        eval_energy_left.push_back( select_EvalEnergy(idbu->methods[coup_left-i-1],idbu->interactions[coup_range-i-1],false) );
//        eval_energy_left.push_back( select_EvalEnergy(idbu->methods[coup_left-i-1],idbu->interactions[coup_left-i-1],false) );
    }
    for (int i=0;i<coup_right;++i) {
        eval_energy_right.push_back( select_EvalEnergy(idbu->methods[coup_right+i+1],idbu->interactions[coup_range+i+1],false) );
//        eval_energy_right.push_back( select_EvalEnergy(idbu->methods[coup_right+i+1],idbu->interactions[coup_right+i+1],false) );
    }
    #else


    #endif


    /*
        Init Ground State
    */
    T0_store   = idbu->Theta0 * M_PI/180;
    T0         = &T0_store;

    R0_store   = getRotMat(T0_store);
    R0_T_store = R0_store.t();
    R0         = &R0_store;
    R0_T       = &R0_T_store;

    if (T0_store(0)==0 && T0_store(1)==0 && T0_store(2)==0) {
        T0_subtract=true;
    }

    local_isotropic = (coupling_range==0 && T0_store(0) == T0_store(1) && check_isotropic_bending());

    // initialize energy pointers as nullptr
    // is properly assigned later
    for (int i=0;i<coup_left;++i) {
        trial_energy_l.push_back(nullptr);
        energy_l.push_back(nullptr);
    }
    for (int i=0;i<coup_right;++i) {
        trial_energy_r.push_back(nullptr);
        energy_r.push_back(nullptr);
    }

    return true;
}


void   BPStep::set_stiffmat() {
    if ( params_diag.size() >= 9 ) {
        params_diag_mat = {{params_diag[0],params_diag[1],params_diag[2]},{params_diag[3],params_diag[4],params_diag[5]},{params_diag[6],params_diag[7],params_diag[8]}};
        params_diag_mat *= STIFFMAT_FACTOR/disc_len;
    }

    params_left_mat = arma::zeros(3,3,coup_left);
    for (unsigned i=0;i<coup_left;i++) {
        if ( params_left[i].size() >= 9 ) {
            params_left_mat.slice(i) = {{params_left[i][0],params_left[i][1],params_left[i][2]},{params_left[i][3],params_left[i][4],params_left[i][5]},{params_left[i][6],params_left[i][7],params_left[i][8]}};
            params_left_mat.slice(i) *= STIFFMAT_FACTOR/disc_len;
        }
    }

    params_right_mat = arma::zeros(3,3,coup_left);
    for (unsigned i=0;i<coup_right;i++) {
        if ( params_right[i].size() >= 9 ) {
            params_right_mat.slice(i) = {{params_right[i][0],params_right[i][1],params_right[i][2]},{params_right[i][3],params_right[i][4],params_right[i][5]},{params_right[i][6],params_right[i][7],params_right[i][8]}};
            params_right_mat.slice(i) *= STIFFMAT_FACTOR/disc_len;
        }
    }
}


////////////////////////////////////////////////////////////////////////////
///////// INIT neighbor BPSteps ////////////////////////////////////////////


bool BPStep::init_neighbors(std::vector<BPStep*>& BPS_all) {
    int len = BPS_all.size();
    int refid;

    for (int nl=0;nl<coup_left;nl++) {
        refid = id-nl-1;
        if (refid<0) {
            if (closed) {
                refid = pmod(refid,len);
            }
            else {
                throw std::domain_error("BPStep::init_neighbors(): Index for neighbour BPStep is out of range. Error in left coupling range.");
                return false;
            }
        }
        BPSl.push_back(BPS_all[refid]);
    }

    for (int nr=0;nr<coup_right;nr++) {
        refid = id+nr+1;
        if (refid>=len) {
            if (closed) {
                refid = pmod(refid,len);
            }
            else {
                throw std::domain_error("BPStep::init_neighbors(): Index for neighbour BPStep is out of range. Error in right coupling range.");
                return false;
            }
        }
        BPSr.push_back(BPS_all[refid]);
    }
    return true;
}

void BPStep::check_neighbors() {
    std::cout << "--------------\n";
    std::cout << "my id " << id << std::endl;
    std::cout << "left neighbor ids" << std::endl;
    for (int nl=0;nl<coup_left;nl++) {
        if (!BPSl[nl]) {
            std::cout << "nullptr!" << std::endl;
        }
        std::cout << "  " << BPSl[nl]->get_id() << std::endl;
    }
    std::cout << "right neighbor ids" << std::endl;
    for (int nr=0;nr<coup_right;nr++) {
        if (!BPSr[nr]) {
            std::cout << "nullptr!" << std::endl;
        }
        std::cout << "  " << BPSr[nr]->get_id() << std::endl;
    }
}

////////////////////////////////////////////////////////////////////////////
/////////// COUPLING MATRIX CONSISTENCY ////////////////////////////////////

/*
Since interaction matrices are assigned to base pair steps, coupling matrices
can potentially be inconsistent which would lead to energetic inconsistencies.
The function BPStep::check_coupling_consistency() checks whether all coupling
matrices are consistent. If wanted consistency can be imposed by using the mean
value of the conflicting matrices. By default the arithmetic mean is calculated
but the the function can be set to calculate the harmonic mean.
In general averaging is not advised and can be avoided by keeping the input files
consistent.
*/

bool BPStep::check_coupling_consistency(bool avg_conflicts, bool harmonic) {

    bool consistent=true;
    for (int nl=0;nl<coup_left;nl++) {

        if ( ! same_EvalEnergy( eval_energy_left[nl] , BPSl[nl]->get_eval_energy_right(nl) )  ) {
            if (avg_conflicts && eval_energy_left[nl]->get_parameter_averaging_allowed()) {
                std::cout << "Warning: Conflict in coupling parameters detected in " << nl+1 << " level left coupling between " << type << " and " << BPSl[nl]->get_type() << std::endl;
                if (!harmonic)  { std::cout << " -> Resolving by taking arithmetic mean." << std::endl; }
                else            { std::cout << " -> Resolving by taking harmonic mean." << std::endl; }

                std::vector<double> avgparams = avg_params(eval_energy_left[nl]->get_params(), BPSl[nl]->get_eval_energy_right(nl)->get_params(),harmonic);
                // set own Ml
                eval_energy_left[nl]->set_params(avgparams);
                // set left neighbors Mr
                BPSl[nl]->get_eval_energy_right(nl)->set_params(avgparams);
            }
            else {
                consistent=false;
            }
        }
    }
    for (int nr=0;nr<coup_right;nr++) {

        if ( ! same_EvalEnergy( eval_energy_right[nr] , BPSr[nr]->get_eval_energy_left(nr) )  ) {
            if (avg_conflicts && eval_energy_right[nr]->get_parameter_averaging_allowed()) {
                std::cout << "Warning: Conflict in coupling parameters detected in " << nr+1 << " level right coupling between " << type << " and " << BPSr[nr]->get_type() << std::endl;
                if (!harmonic)  { std::cout << " -> Resolving by taking arithmetic mean." << std::endl; }
                else            { std::cout << " -> Resolving by taking harmonic mean." << std::endl; }

                std::vector<double> avgparams = avg_params(eval_energy_right[nr]->get_params(), BPSr[nr]->get_eval_energy_left(nr)->get_params(),harmonic);
                // set own Ml
                eval_energy_right[nr]->set_params(avgparams);
                // set left neighbors Mr
                BPSr[nr]->get_eval_energy_left(nr)->set_params(avgparams);
            }
            else {
                consistent=false;
            }
        }
    }
    return consistent;
}


void BPStep::print_M_ptr() {
/*
purely for testing
*/
    for (int nl=coup_left-1;nl>=0;nl--) {
        std::cout << eval_energy_left[nl] << std::endl;
    }
    std::cout << eval_energy_diag << std::endl;
    for (int nr=0;nr<coup_right;nr++) {
        std::cout << eval_energy_right[nr] << std::endl;
    }
}

bool BPStep::same_EvalEnergy(EvalEnergy * eval1,EvalEnergy * eval2) {
    if (eval1->get_method() != eval2->get_method()) {
        return false;
    }
    std::vector<double> params1 = eval1->get_params();
    std::vector<double> params2 = eval2->get_params();
    if (params1.size() != params2.size()) {
        return false;
    }
    for (unsigned i=0;i<params1.size();i++) {
        if (params1[i] != params2[i]) {
            return false;
        }
    }
    return true;
}

std::vector<double> BPStep::avg_params(const std::vector<double> & params1,const std::vector<double> & params2, bool harmonic) {
    std::vector<double> avg_params;
    if (harmonic) {
        for (unsigned i=0;i<params1.size();i++) {
            avg_params.push_back(2*params1[i]*params2[i]/(params1[i]+params2[i]));
        }
    }
    else {
        for (unsigned i=0;i<params1.size();i++) {
            avg_params.push_back(0.5*(params1[i]+params2[i]));
        }
    }
    return avg_params;
}


////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////
////////// ENERGY FUNCTIONS ////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////


/////////////////////////////////////
/////////////// CORE ////////////////


void BPStep::propose_move(const arma::mat& T1, const arma::mat& T2) {
    if (T0_subtract) {
        trial_Theta = ExtractTheta( T1.t() * T2 ) - *T0;
    }
    else {
        trial_Theta = ExtractTheta( *R0_T * T1.t() * T2 );
    }
    trial_pending      = true;
    trial_eval_pending = true;
}

double BPStep::eval_delta_energy() {
/*
This function returns the delta E
*/
    if (!trial_pending || !trial_eval_pending) {
        return 0;
    }

    if (std::abs(trial_Theta(2) - Theta(2)) > M_PI*0.5) {
        /*
            prevent twist discipation by twist exceeding pi.
        */
        return 1e12;
    }

//    trial_energy_0 = arma::dot(trial_Theta, M0 * trial_Theta);
//    trial_energy_0 = eval(trial_Theta);

//    trial_energy_0 = eval_energy_diag->cal_E(trial_Theta);

    trial_energy_0 = eval_energy_diag->cal_beta_energy_diag(trial_Theta);
    double dE = trial_energy_0 - energy_0;
    for (int nl=0;nl<coup_left;nl++) {
        #if PAIR_EVAL_FORWARD == 1
        if ( !( BPSl[nl]->pending() && !BPSl[nl]->eval_pending() ) )
        #else
        if ( !(BPSl[nl]->eval_pending()) )
        #endif
        {
            *trial_energy_l[nl] = eval_energy_left[nl]->cal_beta_energy_offdiag( *(BPSl[nl]->get_trial_Theta()) , trial_Theta);
            /*
            multiplied with two because the other half will not be considered. By assigning this
            value to the pointer it is automatically shared with the coupled BPStep.
            */
            dE += 2*( *trial_energy_l[nl] - *energy_l[nl]);
        }
    }
    for (int nr=0;nr<coup_right;nr++) {
        #if PAIR_EVAL_FORWARD == 1
        if ( !( BPSr[nr]->pending() && !BPSr[nr]->eval_pending() ) )
        #else
        if ( !(BPSr[nr]->eval_pending()) )
        #endif
        {
            *trial_energy_r[nr] = eval_energy_right[nr]->cal_beta_energy_offdiag( trial_Theta , *(BPSr[nr]->get_trial_Theta()) );
            dE += 2*( *trial_energy_r[nr] - *energy_r[nr]);
        }
    }
    trial_eval_pending = false;
    return dE;
}

double BPStep::eval_new_energy() {
/*
This function returns the new energy of the proposed configuration
*/
    if (!trial_pending || !trial_eval_pending) {
        return 0;
    }

//    trial_energy_0 = arma::dot(trial_Theta, M0 * trial_Theta);
//    trial_energy_0 = eval(trial_Theta);

//    trial_energy_0 = eval_energy_diag->cal_E(trial_Theta);

    trial_energy_0 = eval_energy_diag->cal_beta_energy_diag(trial_Theta);
    double dE = trial_energy_0;
    for (int nl=0;nl<coup_left;nl++) {
        #if PAIR_EVAL_FORWARD == 1
        if ( !( BPSl[nl]->pending() && !BPSl[nl]->eval_pending() ) )
        #else
        if ( !(BPSl[nl]->eval_pending()) )
        #endif
        {
            *trial_energy_l[nl] = eval_energy_left[nl]->cal_beta_energy_offdiag( *(BPSl[nl]->get_trial_Theta()) , trial_Theta);
            /*
            multiplied with two because the other half will not be considered. By assigning this
            value to the pointer it is automatically shared with the coupled BPStep.
            */
            dE += 2*( *trial_energy_l[nl]);
        }
    }
    for (int nr=0;nr<coup_right;nr++) {
        #if PAIR_EVAL_FORWARD == 1
        if ( !( BPSr[nr]->pending() && !BPSr[nr]->eval_pending() ) )
        #else
        if ( !(BPSr[nr]->eval_pending()) )
        #endif
        {
            *trial_energy_r[nr] = eval_energy_right[nr]->cal_beta_energy_offdiag( trial_Theta , *(BPSr[nr]->get_trial_Theta()) );
            dE += 2*( *trial_energy_r[nr]);
        }
    }
    trial_eval_pending = false;
    return dE;
}


void BPStep::set_move(bool accepted) {
    trial_pending      = false;
    trial_eval_pending = false;
    if (accepted) {
        Theta    = trial_Theta;
        energy_0 = trial_energy_0;
        for (int nl=0;nl<coup_left;nl++) {
            *energy_l[nl] = *trial_energy_l[nl];
        }
        for (int nr=0;nr<coup_right;nr++) {
            *energy_r[nr] = *trial_energy_r[nr];
        }
    }
}


bool BPStep::eval_pending() {
    return trial_eval_pending;
}

bool BPStep::pending() {
    return trial_pending;
}


/////////////////////////////////////
/////////////////////////////////////
/*
MCStep evaluating elastic energies should
not use these functions. These might be useful
however for extensions that inherit BPStep.
*/

bool BPStep::pending_closed() {
    return (!(trial_eval_pending || trial_pending));
}

void BPStep::close_pending() {
    trial_eval_pending = false;
    trial_pending = false;
}


/////////////////////////////////////
///////////// ACCESSOR //////////////


arma::colvec* BPStep::get_Theta() {
    return &Theta;
}

double BPStep::get_Theta_component(int component) {
    return Theta(component);
}


arma::colvec* BPStep::get_trial_Theta() {
    if (trial_pending) {
        return &trial_Theta;
    }
    else {
        return &Theta;
    }
}


double BPStep::get_energy() {
    double E = energy_0;
    for (int nl=0;nl<coup_left;nl++) {
        E += *energy_l[nl];
    }
    for (int nr=0;nr<coup_right;nr++) {
        E += *energy_r[nr];
    }
    return E;
}

double BPStep::get_trial_energy() {
    if (!trial_pending || trial_eval_pending) {
        return 0;
    }
    double E = trial_energy_0;
    for (int nl=0;nl<coup_left;nl++) {
        E += *trial_energy_l[nl];
    }
    for (int nr=0;nr<coup_right;nr++) {
        E += *trial_energy_r[nr];
    }
    return E;
}


double BPStep::get_energy_single(int displ) {
/*
This function gives access to the individual energies.
*/
    if (displ==0) {
        return energy_0;
    }
    if (displ < 0) {
        if (-displ <= coup_left) {
            return *energy_l[-displ-1];
        }
        else {
            return 0;
        }
    }
    else {
        if (displ <= coup_right) {
            return *energy_r[displ-1];
        }
        else {
            return 0;
        }
    }
}

double BPStep::get_trial_energy_single(int displ) {
/*
This function gives access to the individual energies.
*/
    if (displ==0) {
        return energy_0;
    }
    if (displ < 0) {
        if (-displ <= coup_left) {
            return *energy_l[-displ-1];
        }
        else {
            return 0;
        }
    }
    else {
        if (displ <= coup_right) {
            return *energy_r[displ-1];
        }
        else {
            return 0;
        }
    }
}

//////////////////////////////////////
//////////////////////////////////////
////// Extract Current Energies //////


void    BPStep::energy_extract_select() {
/*
This function selects the energy of the BPStep to be extracted. This function
is important in order to not double count long range interactions.
*/
    energy_extract_selected = true;
}

bool    BPStep::energy_extract_is_selected() {
/*
This function serves as a way for different BPStep to communicate on whether
their energy still has to be extracted. If so the calling BPStep will not
consider this coupling energy.
*/
    return energy_extract_selected;
}

double  BPStep::energy_extract() {
/*
This function returns the current energy of the given BPStep. If inter basepair step
interactions are active and a coupling exists with other BPSteps that are selected
for extraction too, these energy terms will be neglected. This function deactivates
the flag energy_extract_selected so that when BPStep::energy_extract() is called
for the coupled BPSteps they will include the coupling energy with this BPStep.
This method ensures that energy extraction is consistent, i.e. no energies are
considered twice.
*/
    if (!energy_extract_selected) {
        return 0;
    }

    double energy = energy_0;
    for (int nl=0;nl<coup_left;nl++) {
        if ( !BPSl[nl]->energy_extract_is_selected() )
        {
            /*
            multiplied with two because the other half will not be considered.
            */
            energy += 2*( *energy_l[nl] );
        }
    }
    for (int nr=0;nr<coup_right;nr++) {
        if ( !BPSr[nr]->energy_extract_is_selected() )
        {
            /*
            multiplied with two because the other half will not be considered.
            */
            energy += 2*( *energy_r[nr] );
        }
    }
    energy_extract_selected = false;
    return energy;
}

//////////////////////////////////////
//////////////////////////////////////

double BPStep::get_torsional_torque() {
//    cout << -arma::dot(M0_energy->col(2),Theta) << endl;
//    return -arma::dot(M0_energy->col(2),Theta);
}


//////////////////////////////////////
//////////////////////////////////////

bool BPStep::local_isotropic_hamiltonian() {
    return local_isotropic;
}

/////////////////////////////////////
////////// INIT /////////////////////

void BPStep::assign_energy_l_pointers(int id, double* el, double* tel) {
/*
Meant to be called from within another instance of BPStep (init_energies) to wire
the pointers to the energies
*/

    if (energy_l[id]) {
        std::cout << "Warning: BPStep::assign_energy_l_pointers(): energy_l[" << id << "] was already assigned." << std::endl;
    }
    else {
        energy_l[id] = el;
    }
    if (trial_energy_l[id]) {
        std::cout << "Warning: BPStep::assign_energy_l_pointers(): trial_energy_l[" << id << "] was already assigned." << std::endl;
    }
    else {
        trial_energy_l[id] = tel;
    }
}

void BPStep::assign_energy_r_pointers(int id, double* er, double* ter) {
/*
Meant to be called from within another instance of BPStep (init_energies) to wire
the pointers to the energies
*/

    if (energy_r[id]) {
        std::cout << "Warning: BPStep::assign_energy_r_pointers(): energy_r[" << id << "] was already assigned." << std::endl;
    }
    else {
        energy_r[id] = er;
    }
    if (trial_energy_r[id]) {
        std::cout << "Warning: BPStep::assign_energy_r_pointers(): trial_energy_r[" << id << "] was already assigned." << std::endl;
    }
    else {
        trial_energy_r[id] = ter;
    }
}

void BPStep::wire_energies() {

    energy_0        = 0;
    trial_energy_0  = 0;
    trial_pending           = false;
    energy_extract_selected = false;

    for (int nl=0;nl<coup_left;nl++) {
        if (!energy_l[nl]){
            energy_l[nl]        = new double;
            *energy_l[nl]       = 0;
            trial_energy_l[nl]  = new double;
            *trial_energy_l[nl] = 0;
            BPSl[nl]->assign_energy_r_pointers(nl,energy_l[nl],trial_energy_l[nl]);
        }
    }

    for (int nr=0;nr<coup_right;nr++) {
        if (!energy_r[nr]){
            energy_r[nr]        = new double;
            *energy_r[nr]       = 0;
            trial_energy_r[nr]  = new double;
            *trial_energy_r[nr] = 0;
            BPSr[nr]->assign_energy_l_pointers(nr,energy_r[nr],trial_energy_r[nr]);
        }
    }
}




