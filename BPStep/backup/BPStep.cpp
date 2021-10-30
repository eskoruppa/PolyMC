/*
 *  This Class stores the information about a certain base pair step. At this point
 *  energies are assumed to be quadratic in the deformation angles. At a later stage
 *  empirical potentials could be added, also to be stored in this object.
 *
 *  Contains stiffness matrices M0 and various coupling matrices with respect to lower
 *  and upper elements along the chain. All matrices are returned as pointers to arma::mat
 *  objects.
 *      access M0:
 *          - get_M0();
 *          - get_Mc(0);
 *      access lower couplings:
 *          - get_Mc(-dist);
 *      access upper couplings:
 *          - get_Mc(dist);
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
 *  TODO: Some checks in get_Mc() could be omitted if it is ensured that other higher level functions
 *        work properly. However the gain would be minimal
 */

#include "BPStep.h"

BPStep::BPStep(int id_in_chain, const std::string& full_seq, IDB& database, double discretization_length, double Temperature ,bool closed)
: id(id_in_chain), closed(closed), trial_pending(false), disc_len(discretization_length), T(Temperature), T0_subtract(false)
{
    coupling_range = database.get_interaction_range();
    type = find_type(id, full_seq, coupling_range);
    set_params(type, coupling_range, database);
    Theta = arma::zeros(3);
    zeromat = arma::zeros(3,3);
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

arma::mat * BPStep::get_M0() {
    return M0;
}

arma::mat * BPStep::get_Ml(int nl) {
    if (nl < coup_left) {
        return Ml[nl];
    }
    else {
        return &zeromat;
    }
}

arma::mat * BPStep::get_Mr(int nr) {
    if (nr < coup_right) {
        return Mr[nr];
    }
    else {
        return &zeromat;
    }
}


arma::mat * BPStep::get_Mc(int displ) {
    if (displ == 0) {
        return M0;
    }
    if (displ < 0) {
        if (-displ <= coup_left) {
            return Ml[-displ-1];
            //  return (Ml+(-displ-1));
        }
        else {
            throw std::invalid_argument( "BPStep::get_MC - Basepair Step coupling out of range" );
            return &zeromat;
        }
    }
    else {
        if (displ <= coup_right) {
            return Mr[displ-1];
            //  return (Mu+(displ-1));
        }
        else {
            throw std::invalid_argument( "BPStep::get_MC - Basepair Step coupling out of range" );
            return &zeromat;
        }
    }
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


void BPStep::change_T(double Tnew) {
    double resc_fac = T/Tnew;
    M0_store = resc_fac*M0_store;
    for (int i=0;i<coup_left;++i) {
        Ml_store[i] = resc_fac*Ml_store[i];
    }
    for (int i=0;i<coup_right;++i) {
        Mr_store[i] = resc_fac*Mr_store[i];
    }
    T = Tnew;
}

void BPStep::set_T0_subtract(bool subtract) {
    T0_subtract = subtract;
}


////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////
///////////// INIT METHODS /////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////

std::string BPStep::find_type(int bps, const std::string& full_seq, int coup_range) {
    int      num_bp = int(full_seq.length());
    int      num_bps;
    std::string   assign_type;
    int      shift;
    int      seglen = (1+coup_range)*2;

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
                if (print_map) {cout << bps << " " << assign_type << "   <--" << shift << endl;}
            }
            // close to right terminus
            if (int(bps+coup_range) >= num_bps) {
                shift = bps+coup_range+1-num_bps;
                assign_type = full_seq.substr(bps-coup_range,seglen-shift);
                for (int i=0;i<shift;++i) {
                    assign_type=assign_type+"x";
                }
                if (print_map) {cout << bps << " " << assign_type << "   <--" << shift << endl;}
            }
        }
        // central
        else {
            assign_type = full_seq.substr(bps-coup_range,seglen);
            if (print_map) {cout << bps << " " << assign_type << endl;}
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
                if (print_map) {cout << bps << " " << assign_type << "   <--" << shift << endl;}
            }
            // close to right terminus
            if (int(bps+coup_range) >= (num_bps-1)) {
                shift = bps+coup_range+2-num_bps;
                assign_type = full_seq.substr(bps-coup_range,seglen-shift);
                assign_type = assign_type+full_seq.substr(0,shift);
                if (print_map) {cout << bps << " " << assign_type << "   <--" << shift << endl;}
            }
        }
        // central
        else {
            assign_type = full_seq.substr(bps-coup_range,seglen);
            if (print_map) {cout << bps << " " << assign_type << endl;}
        }
    }
    return assign_type;
}

////////////////////////////////////////////////////////////////////////////
///////////////// ASSIGN INTERACTION PARAMETERS ////////////////////////////

bool BPStep::set_params(const string& type, int coup_range, IDB& data) {

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
    IDBUnit* idbu = new IDBUnit;

    if (terminus){
        if (!data.IDBU_get_terminus(*idbu, type)) {
            throw std::domain_error("BPStep::set_params(): terminus could not be found");
            return false;
        }
    }
    else {
        idbu = data.IDBU_get_unit(type);
    }

    /*
    This factor is here already being included in the stiffness matrices to avoid
    having to consider it at every evaluations of the stiffness matrix.
    For the neighbor coupling matrices an additional factor of 0.5 is introduced
    to account for the fact that those contributions appear both in this instance of
    BPStep as well as in the paired neighbor.
    */

    double E_fac = 0.5/disc_len * T_ref/T;

    M0_store = E_fac*idbu->Ms.slice(coup_range);
    M0       = &M0_store;

    M0_energy_store = idbu->Ms.slice(coup_range)*4.114/disc_len;
    M0_energy       = &M0_energy_store;

    T0_store = idbu->Theta0 * M_PI/180;
    T0       = &T0_store;

    R0_store   = getRotMat(T0_store);
    R0_T_store = R0_store.t();
    R0   = &R0_store;
    R0_T = &R0_T_store;

    for (int i=0;i<coup_left;++i) {
        Ml_store.push_back(idbu->Ms.slice(coup_left-i-1));
    }
    for (int i=0;i<coup_left;++i) {
        Ml_store[i] = 0.5*E_fac*Ml_store[i];
        Ml.push_back(&Ml_store[i]);
    }

    for (int i=0;i<coup_right;++i) {
        Mr_store.push_back(idbu->Ms.slice(coup_right+i+1));
    }
    for (int i=0;i<coup_right;++i) {
        Mr_store[i] = 0.5*E_fac*Mr_store[i];
        Mr.push_back(&Mr_store[i]);
    }

    // initialize energy pointers as nullptr
    for (int i=0;i<coup_left;++i) {
        trial_energy_l.push_back(nullptr);
        energy_l.push_back(nullptr);
    }
    for (int i=0;i<coup_right;++i) {
        trial_energy_r.push_back(nullptr);
        energy_r.push_back(nullptr);
    }

    if (terminus) {
        delete idbu;
    }
    return true;
}

////////////////////////////////////////////////////////////////////////////
///////// INIT neighbor BPSteps ////////////////////////////////////////////


bool BPStep::init_neighbors(vector<BPStep*>& BPS_all) {
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
    cout << "--------------\n";
    cout << "my id " << id << endl;
    cout << "left neighbor ids" << endl;
    for (int nl=0;nl<coup_left;nl++) {
        if (!BPSl[nl]) {
            cout << "nullptr!" << endl;
        }
        cout << "  " << BPSl[nl]->get_id() << endl;
    }
    cout << "right neighbor ids" << endl;
    for (int nr=0;nr<coup_right;nr++) {
        if (!BPSr[nr]) {
            cout << "nullptr!" << endl;
        }
        cout << "  " << BPSr[nr]->get_id() << endl;
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
        if (!same_mat(*Ml[nl], *BPSl[nl]->get_Mr(nl) ) ) {
            if (avg_conflicts) {
                cout << "Warning: Conflict in Coupling Matrix detected in " << nl+1 << " level left coupling between " << type << " and " << BPSl[nl]->get_type() << endl;
                if (!harmonic)  { cout << " -> Resolving by taking arithmetic mean." << endl; }
                else            { cout << " -> Resolving by taking harmonic mean." << endl; }
                arma::mat avgM = avg_mat(*Ml[nl],*BPSl[nl]->get_Mr(nl),harmonic);
                // set own Ml
                set_Ml(nl,avgM);
                // set left neighbors Mr
                BPSl[nl]->set_Mr(nl,avgM);
            }
            consistent=false;
        }
    }
    for (int nr=0;nr<coup_right;nr++) {
        if (!same_mat(*Mr[nr], *BPSr[nr]->get_Ml(nr) ) ) {
            if (avg_conflicts) {
                cout << "Warning: Conflict in Coupling Matrix detected in " << nr+1 << " level right coupling between " << type << " and " << BPSr[nr]->get_type() << endl;
                if (!harmonic)  { cout << " -> Resolving by taking arithmetic mean." << endl; }
                else            { cout << " -> Resolving by taking harmonic mean." << endl; }
                arma::mat avgM = avg_mat(*Mr[nr],*BPSr[nr]->get_Ml(nr),harmonic);
                // set own Ml
                set_Mr(nr,avgM);
                // set right neighbors Ml
                BPSr[nr]->set_Ml(nr,avgM);
            }
            consistent=false;
        }
    }
    return consistent;
}


void BPStep::print_M_ptr() {
/*
purely for testing
*/
    for (int nl=coup_left-1;nl>=0;nl--) {
        cout << Ml[nl] << endl;
    }
    cout << M0 << endl;
    for (int nr=0;nr<coup_right;nr++) {
        cout << Mr[nr] << endl;
    }
}


void BPStep::set_Ml(int nl, arma::mat & M) {
    if (nl < coup_left) {
        Ml_store[nl] = M;
    }
}

void BPStep::set_Mr(int nr, arma::mat & M) {
    if (nr < coup_right) {
        Mr_store[nr] = M;
    }
}

bool BPStep::same_mat(arma::mat & M1, arma::mat & M2) {
    double eps = 1e-10;
    if (arma::accu(abs(M1-M2))<eps) {
        return true;
    }
    return false;
}

arma::mat BPStep::avg_mat(arma::mat & M1, arma::mat & M2, bool harmonic) {
    unsigned dim = M1.n_rows;
    arma::mat avgM = arma::zeros(dim,dim);
    if (harmonic) {
        for (unsigned c=0;c<dim;c++) {
            for (unsigned r=0;r<dim;r++) {
                avgM(r,c) = 2*M1(r,c)*M2(r,c)/( M1(r,c)+M2(r,c) );
            }
        }
    }
    else {
        avgM = 0.5*(M1+M2);
    }
    return avgM;
}


////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////
////////// ENERGY FUNCTIONS ////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////


/////////////////////////////////////
/////////////// CORE ////////////////

//void BPStep::propose_move(arma::colvec& Theta) {
////    arma::colvec T = *T0;
////    arma::mat Rz = Rotz( T(2));
////    arma::mat R0 = getRotMat( *T0 );
////    arma::mat R = getRotMat(Theta);
////    trial_Theta = ExtractTheta(R0.t()*R);
//    trial_Theta        = Theta-*T0;
//    trial_pending      = true;
//    trial_eval_pending = true;
//}

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

    trial_energy_0 = arma::dot(trial_Theta, *M0 * trial_Theta);
    double dE = trial_energy_0 - energy_0;
    for (int nl=0;nl<coup_left;nl++) {
        #if PAIR_EVAL_FORWARD == 1
        if ( !( BPSl[nl]->pending() && !BPSl[nl]->eval_pending() ) )
        #else
        if ( !(BPSl[nl]->eval_pending()) )
        #endif
        {
            *trial_energy_l[nl] = arma::dot( *(BPSl[nl]->get_trial_Theta()) , *(Ml[nl]) * trial_Theta );
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
            *trial_energy_r[nr] = arma::dot( trial_Theta, *(Mr[nr]) * *(BPSr[nr]->get_trial_Theta()) );
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

    trial_energy_0 = arma::dot(trial_Theta, *M0 * trial_Theta);
    double dE = trial_energy_0;
    for (int nl=0;nl<coup_left;nl++) {
        #if PAIR_EVAL_FORWARD == 1
        if ( !( BPSl[nl]->pending() && !BPSl[nl]->eval_pending() ) )
        #else
        if ( !(BPSl[nl]->eval_pending()) )
        #endif
        {
            *trial_energy_l[nl] = arma::dot( *(BPSl[nl]->get_trial_Theta()) , *(Ml[nl]) * trial_Theta );
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
            *trial_energy_r[nr] = arma::dot( trial_Theta, *(Mr[nr]) * *(BPSr[nr]->get_trial_Theta()) );
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
    return -arma::dot(M0_energy->col(2),Theta);
}


//////////////////////////////////////
//////////////////////////////////////

bool BPStep::local_isotropic_hamiltonian() {
    if ( coupling_range > 0             ||
         M0_store(0,0) != M0_store(1,1) ||
         M0_store(0,1) != 0             ||
         M0_store(0,2) != 0             ||
         M0_store(1,2) != 0 )
    {
        return false;
    }
    return true;
}




/////////////////////////////////////
////////// INIT /////////////////////

void BPStep::assign_energy_l_pointers(int id, double* el, double* tel) {
/*
Meant to be called from within another instance of BPStep (init_energies) to wire
the pointers to the energies
*/

    if (energy_l[id]) {
        cout << "Warning: BPStep::assign_energy_l_pointers(): energy_l[" << id << "] was already assigned." << endl;
    }
    else {
        energy_l[id] = el;
    }
    if (trial_energy_l[id]) {
        cout << "Warning: BPStep::assign_energy_l_pointers(): trial_energy_l[" << id << "] was already assigned." << endl;
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
        cout << "Warning: BPStep::assign_energy_r_pointers(): energy_r[" << id << "] was already assigned." << endl;
    }
    else {
        energy_r[id] = er;
    }
    if (trial_energy_r[id]) {
        cout << "Warning: BPStep::assign_energy_r_pointers(): trial_energy_r[" << id << "] was already assigned." << endl;
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



