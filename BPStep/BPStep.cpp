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
    set_params(type, coupling_range, database);
    Theta = arma::zeros(3);

    #ifdef INCLUDE_GENERALIZED_FORCE
    generalized_force_active    = false;
    generalized_force           = arma::zeros(3);

    generalied_force_energy          = 0;
    generalied_force_trial_energy    = 0;
    #endif

    #ifdef BPS_USE_EVALENERGY
    cov = eval_energy_diag->get_cov();
    #else
    cov = cal_cov();
    #endif

    kT = ref_kT*temp/ref_temp;
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
    return cov;
}

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

    #ifdef BPS_USE_EVALENERGY
    eval_energy_diag->set_temp(new_temp);
    for (int i=0;i<coup_left;++i) {
        eval_energy_left[i]->set_temp(new_temp);
    }
    for (int i=0;i<coup_right;++i) {
        eval_energy_right[i]->set_temp(new_temp);
    }
    temp = new_temp;
    cov = eval_energy_diag->get_cov();

    #else
    params_diag_mat = params_diag_mat*temp/new_temp;
    for (unsigned i=0;i<coup_left;i++) {
        params_left_mat.slice(i) *= new_temp/temp;
    }
    for (unsigned i=0;i<coup_right;i++) {
        params_right_mat.slice(i) *= new_temp/temp;
    }

    params_diag      = set_params_temp(new_temp,coup_diag_method, params_diag);
    params_diag_mat  = vec2mat(params_diag);
    for (unsigned i=0;i<coup_left;i++) {
        params_left[i]           = set_params_temp(new_temp,coup_left_method[i], params_left[i]);
        params_left_mat.slice(i) = vec2mat(params_left[i]);
    }
    for (unsigned i=0;i<coup_right;i++) {
        params_right[i] = set_params_temp(new_temp,coup_right_method[i], params_right[i]);
        params_right_mat.slice(i) = vec2mat(params_right[i]);
    }
    temp = new_temp;
    kT   = ref_kT*temp/ref_temp;
    cov  = cal_cov();
    #endif
}

void BPStep::set_T0_subtract(bool subtract) {
    T0_subtract = subtract;
    if (T0_store(0)==0 && T0_store(1)==0 && T0_store(2)==0) {
        T0_subtract=true;
    }
}

void BPStep::deactivate_twist_energy() {
    #ifdef BPS_USE_EVALENERGY
    eval_energy_diag->deactivate_twist_energy();
    for (int i=0;i<coup_left;++i) {
        eval_energy_left[i]->deactivate_twist_energy();
    }
    for (int i=0;i<coup_right;++i) {
        eval_energy_right[i]->deactivate_twist_energy();
    }
    #else
    std::throw( "Error: The option to de-/reactivate the twist energy is currently only supported by the EE_StiffMat object");
    #endif
    twist_energy_active = false;
}
void BPStep::reactivate_twist_energy() {
    #ifdef BPS_USE_EVALENERGY
    eval_energy_diag->reactivate_twist_energy();
    for (int i=0;i<coup_left;++i) {
        eval_energy_left[i]->reactivate_twist_energy();
    }
    for (int i=0;i<coup_right;++i) {
        eval_energy_right[i]->reactivate_twist_energy();
    }
    #else
    std::throw( "Error: The option to de-/reactivate the twist energy is currently only supported by the EE_StiffMat object");
    #endif
    twist_energy_active = true;
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


    #ifdef BPS_USE_EVALENERGY
    /*
        Init energy evaluation classes derived from EvalEnergy
    */

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
    /*
        Init parameters for internal energy calculation
    */
    /*
        This seems weird... set_params???
    */

    coup_diag_method = get_method_id(idbu->methods[coup_range]);
    params_diag      = set_params(coup_diag_method, idbu->interactions[coup_range]);
    params_diag_mat  = vec2mat(params_diag);

    params_left_mat = arma::cube(3,3,coup_range);
    for (int i=0;i<coup_left;++i) {
        coup_left_method.push_back( get_method_id(idbu->methods[coup_left-i-1]) );
        params_left     .push_back( set_params(coup_left_method[i] , idbu->interactions[coup_range-i-1]) );
        params_left_mat.slice(i)  = vec2mat(params_left[i]);
    }
    params_right_mat = arma::cube(3,3,coup_range);
    for (int i=0;i<coup_right;++i) {
        coup_right_method.push_back( get_method_id(idbu->methods[coup_right+i+1]) );
        params_right     .push_back( set_params(coup_right_method[i], idbu->interactions[coup_range+i+1] ) );
        params_right_mat.slice(i) = vec2mat(params_right[i]);
    }
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

    #ifdef BPS_USE_EVALENERGY
    local_isotropic = (coupling_range==0 && T0_store(0) == 0 && T0_store(1) == 0 && eval_energy_diag->isotropic_bending());
    #else
    local_isotropic = (coupling_range==0 && T0_store(0) == 0 && T0_store(1) == 0 && check_isotropic_bending());
    #endif


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

arma::mat BPStep::vec2mat(std::vector<double> params) {
    if (params.size() >= 9) {
        return {{params[0],params[1],params[2]},{params[3],params[4],params[5]},{params[6],params[7],params[8]}};
    }
    return arma::zeros(3,3);
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
    #ifdef BPS_USE_EVALENERGY
    for (int nl=0;nl<coup_left;nl++) {

        if ( ! same_EvalEnergy( eval_energy_left[nl] , BPSl[nl]->get_eval_energy_right(nl) )  ) {

            std::cout << "Warning: Conflict in coupling parameters detected in " << nl+1 << " level left coupling between " << type << " and " << BPSl[nl]->get_type() << std::endl;
            if (avg_conflicts && eval_energy_left[nl]->get_parameter_averaging_allowed()) {

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

            std::cout << "Warning: Conflict in coupling parameters detected in " << nr+1 << " level right coupling between " << type << " and " << BPSr[nr]->get_type() << std::endl;
            if (avg_conflicts && eval_energy_right[nr]->get_parameter_averaging_allowed()) {

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
    #else

    for (int nl=0;nl<coup_left;nl++) {
        if (! equal_params( params_left[nl] ,  BPSl[nl]->params_right[nl] ) ) {

            std::cout << "Warning: Conflict in coupling parameters detected in " << nl+1 << " level left coupling between " << type << " and " << BPSl[nl]->get_type() << std::endl;
            if (avg_conflicts && allow_params_averaging(coup_left_method[nl])) {

                if (!harmonic)  { std::cout << " -> Resolving by taking arithmetic mean." << std::endl; }
                else            { std::cout << " -> Resolving by taking harmonic mean." << std::endl; }

                std::vector<double> avgparams = avg_params(params_left[nl], BPSl[nl]->params_right[nl],harmonic);
                // set own
                params_left[nl] = avgparams;
                // set left neighbors
                BPSl[nl]->params_right[nl] = avgparams;
            }
            else {
                consistent=false;
            }
        }
    }
    for (int nr=0;nr<coup_right;nr++) {

        if (! equal_params( params_right[nr] ,  BPSr[nr]->params_left[nr] ) ) {

            std::cout << "Warning: Conflict in coupling parameters detected in " << nr+1 << " level right coupling between " << type << " and " << BPSr[nr]->get_type() << std::endl;
            if (avg_conflicts && allow_params_averaging(coup_right_method[nr])) {

                if (!harmonic)  { std::cout << " -> Resolving by taking arithmetic mean." << std::endl; }
                else            { std::cout << " -> Resolving by taking harmonic mean." << std::endl; }

                std::vector<double> avgparams = avg_params(params_right[nr], BPSr[nr]->params_left[nr],harmonic);
                // set own
                params_right[nr] = avgparams;
                // set left neighbors
                BPSr[nr]->params_left[nr] = avgparams;
            }
            else {
                consistent=false;
            }
        }
    }

    #endif
    return consistent;
}


void BPStep::print_M_ptr() {
/*
purely for testing
*/
    #ifdef BPS_USE_EVALENERGY
    for (int nl=coup_left-1;nl>=0;nl--) {
        std::cout << eval_energy_left[nl] << std::endl;
    }
    std::cout << eval_energy_diag << std::endl;
    for (int nr=0;nr<coup_right;nr++) {
        std::cout << eval_energy_right[nr] << std::endl;
    }
    #endif
}

bool BPStep::equal_params(const std::vector<double> & params1,const std::vector<double> & params2) {
    if (params1.size() != params2.size()) {
        return false;
    }
    for (unsigned i=0;i<params1.size();i++) {
        if (params1[i] != params2[i]) return false;
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

    if (twist_energy_active && std::abs(trial_Theta(2) - Theta(2)) > M_PI*0.5) {
        /*
            prevent twist discipation by twist exceeding pi.
        */
        return 1e12;
    }

    #ifdef BPS_USE_EVALENERGY

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

    #else
    trial_energy_0 = arma::dot(trial_Theta, params_diag_mat * trial_Theta);

//    trial_energy_0 = cal_beta_energy_diag(trial_Theta);
    double dE = trial_energy_0 - energy_0;
    for (int nl=0;nl<coup_left;nl++) {
        #if PAIR_EVAL_FORWARD == 1
        if ( !( BPSl[nl]->pending() && !BPSl[nl]->eval_pending() ) )
        #else
        if ( !(BPSl[nl]->eval_pending()) )
        #endif
        {
            *trial_energy_l[nl] = cal_beta_energy_left(*(BPSl[nl]->get_trial_Theta()),trial_Theta,nl);
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
            *trial_energy_r[nr] = cal_beta_energy_right(trial_Theta,*(BPSr[nr]->get_trial_Theta()),nr);
            dE += 2*( *trial_energy_r[nr] - *energy_r[nr]);
        }
    }

    #endif

    #ifdef INCLUDE_GENERALIZED_FORCE
    if (generalized_force_active) {
        generalied_force_trial_energy = eval_generalized_force(trial_Theta);
        dE += generalied_force_trial_energy - generalied_force_energy;
    }
    #endif

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

    #ifdef BPS_USE_EVALENERGY

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

    #else
    trial_energy_0 = arma::dot(trial_Theta, params_diag_mat * trial_Theta);

//    trial_energy_0 = cal_beta_energy_diag(trial_Theta);
    double dE = trial_energy_0;
    for (int nl=0;nl<coup_left;nl++) {
        #if PAIR_EVAL_FORWARD == 1
        if ( !( BPSl[nl]->pending() && !BPSl[nl]->eval_pending() ) )
        #else
        if ( !(BPSl[nl]->eval_pending()) )
        #endif
        {
            *trial_energy_l[nl] = cal_beta_energy_left(*(BPSl[nl]->get_trial_Theta()),trial_Theta,nl);
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
            *trial_energy_r[nr] = cal_beta_energy_right(trial_Theta,*(BPSr[nr]->get_trial_Theta()),nr);
            dE += 2*( *trial_energy_r[nr]);
        }
    }

    #endif

    #ifdef INCLUDE_GENERALIZED_FORCE
    if (generalized_force_active) {
        generalied_force_trial_energy = eval_generalized_force(trial_Theta);
        dE += generalied_force_trial_energy;
    }
    #endif

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

        #ifdef INCLUDE_GENERALIZED_FORCE
        if (generalized_force_active) {
            generalied_force_energy = generalied_force_trial_energy;
        }
        #endif
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

arma::colvec BPStep::get_FullTheta() {
    return Theta + (*T0);
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
    std::cout << "ERROR: function BPStep::get_torsional_torque() not yet implemented!" << std::endl;
    std::exit(0);
    return 0;
}


//////////////////////////////////////
//////////////////////////////////////

bool BPStep::local_isotropic_hamiltonian() {
    return local_isotropic;
}

std::vector<double> BPStep::get_status_diag() {
    #ifdef BPS_USE_EVALENERGY
    return eval_energy_diag->get_status_diag(Theta);
    #else
    std::cout << "Warning: BPStep::get_status_diag() only works if BPS_USE_EVALENERGY is selected!" std::endl;
    std::cout << " -> returning empty status vector" std::endl;
    return {};
    #endif
}

std::string * BPStep::get_method_diag() {
    #ifdef BPS_USE_EVALENERGY
    return eval_energy_diag->get_method();
    #else
    std::cout << "Warning: BPStep::get_method_diag() only works if BPS_USE_EVALENERGY is selected!" std::endl;
    std::cout << " -> returning NULL pointer" std::endl;
    return NULL;
    #endif
}

EvalEnergy * BPStep::get_EvalEnergy_diag() {
    return eval_energy_diag;
}


/////////////////////////////////////////////
arma::colvec BPStep::Triads2Theta       (const arma::mat T1,const arma::mat T2) {
    if (T0_subtract) {
        return ExtractTheta( T1.t() * T2 ) - *T0;
    }
    else {
        return ExtractTheta( *R0_T * T1.t() * T2 );
    }
}
arma::colvec BPStep::Triads2FullTheta   (const arma::mat T1,const arma::mat T2){
    return ExtractTheta( T1.t() * T2 );
}
arma::mat    BPStep::Theta2Rotmat       (const arma::colvec Theta) {
    if (T0_subtract) {
        return getRotMat(Theta + *T0);
    }
    else {
        return *R0 * getRotMat(Theta);
    }
}
arma::mat    BPStep::FullTheta2Rotmat   (const arma::colvec FullTheta) {
    return getRotMat(FullTheta);
}


//////////////////////////////////////
//////////////////////////////////////
/////// GENERALIZED FORCE ////////////

#ifdef INCLUDE_GENERALIZED_FORCE
void BPStep::set_generalized_force(const arma::colvec & gen_force) {

    if (arma::norm(gen_force) > 1e-10) {
        generalized_force_active = true;
        generalized_force        = gen_force;
    }
    else {
        deactivate_generalized_force();
    }
}

void BPStep::deactivate_generalized_force(){
    generalized_force_active = false;
    generalized_force        = arma::zeros(3);
}

double BPStep::eval_generalized_force(const arma::colvec & Theta){
    return -arma::dot(generalized_force,Theta)/kT;
}

double BPStep::eval_change_generalized_force(const arma::colvec & Theta){
    return arma::dot(generalized_force,Theta)/kT - generalied_force_energy;
}
#endif


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


///////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////
////////////////////// EVAL ENERGY ////////////////////////////////////////////////////////

#ifdef BPS_USE_EVALENERGY
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
#endif


std::vector<double> BPStep::get_eval_energy_params(int id) {
    #ifdef BPS_USE_EVALENERGY

    if (id > 0) {
        int nr = id-1;
        if (nr >= coup_right) {
            return {};
        }
        else {
            return eval_energy_right[nr]->get_params();
        }
    }
    if (id < 0) {
        int nl = -id-1;
        if (nl >= coup_left) {
            return {};
        }
        else {
            return eval_energy_left[nl]->get_params();
        }
    }
    return eval_energy_diag->get_params();

    #else
    std::cout << "Error: BPStep::get_eval_energy_params() is only defined for EVALENERGY." << std::endl;
    std::exit(0);
    #endif
}

