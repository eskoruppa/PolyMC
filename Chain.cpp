#include "Chain.h"

/*


*/

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////// CONSTRUCTOR ///////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

Chain::Chain(std::string interaction_file)
: data(interaction_file), conf_initialized(false) {
    set_interaction_range(data.get_interaction_range());
    disc_len = data.get_disc_len();
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////// DESTRUCTOR ///////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

Chain::~Chain() {
}


////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////// ACCESSORS ////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

std::vector<BPStep*>* Chain::get_BPS() {
    return &BPS;
}
arma::mat*  Chain::get_bp_pos() {
    return &bp_pos;

}
arma::cube* Chain::get_triads() {
    return &triads;
}

double Chain::get_disc_len() {
    return disc_len;
}

double Chain::get_contour_len() {
    return contour_len;
}

int Chain::get_num_bp() {
    return num_bp;
}

int Chain::get_num_bps() {
    return num_bps;
}

bool Chain::local_isotropic_hamiltonian() {
    for (unsigned i=0;i<num_bps;i++) {
        if (!BPS[i]->local_isotropic_hamiltonian()) {
            return false;
        }
    }
    return true;
}

double   Chain::get_kT() {
    return kT;
}

double   Chain::get_beta() {
    return beta;
}

double   Chain::get_T() {
    return T;
}
double   Chain::get_T_ref() {
    return T_ref;
}

unsigned Chain::get_interaction_range() {
    return interaction_range;
}

void Chain::cal_avg_stiff(bool inittrue) {
    if (conf_initialized || inittrue) {
        avg_cov = arma::zeros(3,3);
        for (unsigned i=0;i<num_bps;i++) {
            avg_cov += BPS[i]->get_cov();
        }
        avg_cov   = avg_cov/num_bps;
        avg_stiff = arma::inv(avg_cov)*disc_len;
        avg_chol  = arma::chol(avg_cov,"lower");
    }
}

arma::mat Chain::get_avg_stiff() {
    if (conf_initialized) {
        return avg_stiff;
    }
    else {
        return arma::zeros(3,3);
    }
}

arma::mat Chain::get_avg_cov() {
/*
    Returns the average covariance matrix of the on site coupling matrices. This
    matrix is dimensionless, i.e. it is the inverse average stiffness matrix times the
    discretization length.
*/
    if (conf_initialized) {
        return avg_cov;
    }
    else {
        return arma::zeros(3,3);
    }
}

arma::mat Chain::get_avg_chol() {
/*
    Returns the lower triangular component of the choleski decomposition of the
    average covariance matrix.
*/
    if (conf_initialized) {
        return avg_chol;
    }
    else {
        return arma::zeros(3,3);
    }
}


bool     Chain::force_active(){
    return force_constrained;
}
double   Chain::get_force() {
    return force;
}
arma::colvec Chain::get_force_dir() {
    return force_dir;
}
arma::colvec Chain::get_beta_force_vec() {
    return beta_force_vec;
}
bool     Chain::fixed_termini() {
    return termini_fixed;
}
bool     Chain::fixed_termini_radial() {
    return termini_fix_radial;
}
bool     Chain::fixed_first_orientation() {
    return first_fixed_orientation;
}
bool     Chain::fixed_last_orientation() {
    return last_fixed_orientation;
}
bool     Chain::topology_closed() {
    return closed_topology;
}
std::string Chain::get_config_type() {
    return config_type;
}
std::string Chain::get_sequence() {
    return seq;
}



////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////// SETTINGS MUTATORS ///////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/*
All setting are only change in the mutator methods. Other
member functions that require the change of some of the
settings call these methods to keep the setting consistent
*/

void Chain::set_helical_repeat_length(double hel_rep){
    helical_repeat_length       = hel_rep;
    avg_intrinsic_twist_density = 2*M_PI/helical_repeat_length;
}

void Chain::set_intrinsic_twist_density(double twist_density){
    avg_intrinsic_twist_density = twist_density;
    helical_repeat_length       = 2*M_PI/avg_intrinsic_twist_density;
}

void Chain::set_T(double temp) {
    T        = temp;
    kT       = kT_ref*temp/T_ref;
    beta     = 1./kT;
    invT_fac = T_ref/T;

    change_torque (torque);
    set_force  (force,force_dir);

    if (conf_initialized) {
        for (unsigned bps=0;bps<num_bps;bps++) {
            BPS[bps]->change_T(temp);
        }
    }
    cal_avg_stiff();
    std::cout << "Temperature set to " << T << "K" << std::endl;
}

void Chain::set_T0_subtract(bool subtract) {
    T0_subtract = subtract;
    if (conf_initialized) {
        for (unsigned bps=0;bps<num_bps;bps++) {
            BPS[bps]->set_T0_subtract(T0_subtract);
        }
    }
}

void Chain::set_interaction_range(unsigned ir) {
    if (ir > data.get_interaction_range()) {
        std::cout << "Warning: Interaction range change ignored. The provided interaction database only contains" << std::endl;
        std::cout << "interactions up to a range of " << data.get_interaction_range() << " (" << ir << " given)" << std::endl;
    }
    else {
        interaction_range = ir;
    }
}


///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


void Chain::set_force(double f, const arma::colvec& dir) {
    if (f==0) { force_constrained = false; }
    else      { force_constrained = true; }
    force             = f;
    if (force==0) force_dir = {0,0,1};
    else          force_dir = dir/arma::norm(dir);
    beta_force_vec    = f*force_dir/kT;
//    cout << "Force set to " << f << " pN" << std::endl;
}

void Chain::fix_termini_position(bool fix_radial, bool fix) {
    if (closed_topology) {
        termini_fixed      = false;
        termini_fix_radial = false;
    }
    else {
    termini_fixed       = fix;
    termini_fix_radial  = fix_radial;
    }
}

void Chain::fix_termini_orientation(bool fix) {
    if (closed_topology) {
        last_fixed_orientation  = false;
        first_fixed_orientation = false;
    }
    else {
        last_fixed_orientation  = fix;
        first_fixed_orientation = fix;
    }
}


////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////// PRIVATE MUTATORS ////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void Chain::set_closed_topology(bool closed) {
    closed_topology = closed;

    if (closed) {
        if (last_fixed_orientation) {
            std::cout << "Warning: Fix last bead orientation overwritten by closed topology." << std::endl;
        }
        if (first_fixed_orientation) {
            std::cout << "Warning: Fix first bead orientation overwritten by closed topology." << std::endl;
        }
        fix_termini_orientation(false);
        if (termini_fixed) {
            std::cout << "Warning: Fix termini overwritten by closed topology." << std::endl;
        }
        fix_termini_position(false,false);
        if (force_constrained) {
            std::cout << "Warning: Force deactivated due to closed topology." << std::endl;
        }
        set_force(0);
        if (link_torsional_trap || link_const_torque) {
            std::cout << "Warning: Torque deactivated due to closed topology." << std::endl;
        }
        change_torque(0);
    }
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////// CONFIGURATION INITIALIZATION //////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

bool Chain::gen_linear(unsigned number_bp, const std::string& sequence , double supercoiling_density, const arma::colvec direction)
{
    std::cout << "Initializing linear chain configuration" << std::endl;
    if (conf_initialized) {
        std::cout << "Warning: Trying to generate linear configuration, but a configuration is already initialized" << std::endl;
        return false;
    }

    num_bp      = number_bp;
    num_bps     = number_bp-1;
    disc_len    = data.get_disc_len();
    contour_len = num_bps*disc_len;
    set_seq(sequence,num_bp);
    set_closed_topology(false);
    config_type = "linear";

    // Initialize BPStep
    init_BPS();
    dLK = sigma2dLk(supercoiling_density);

    // Generate Configuration
    bp_pos = arma::zeros(3,num_bp);
    triads = arma::zeros(3,3,num_bp);

    double exc_tw_pbp = dLK/num_bps*2*M_PI;
    arma::mat RzN     = arma::eye(3,3);
    arma::mat Rz;
    arma::colvec T0;

    arma::colvec initdir = direction/norm(direction);
    arma::colvec refdir  = {0,0,1};

    arma::mat rot_mat;
    if (arma::norm(refdir-initdir)>1e-10) {
        arma::colvec axis = arma::cross(refdir,initdir);
        double theta      = asin(arma::norm(axis));
        if (arma::dot(refdir,initdir)<0) {
            theta = M_PI-theta;
        }
        axis    = theta*axis/arma::norm(axis);
        rot_mat = getRotMat(axis);
    }
    else {
        rot_mat = arma::eye(3,3);
    }

	for (unsigned bp=0; bp<num_bp;bp++) {
		bp_pos.col(bp)   = initdir*disc_len*bp;
		triads.slice(bp) = rot_mat*RzN*arma::eye(3,3);

		double diff = arma::norm(initdir-triads.slice(bp).col(2));
        if (diff > 1e-10) {
            std::cout << "INIT ERROR: WRONG ROTATION!" << std::endl;
            std::cout << diff << std::endl;
            std::cout << arma::norm(initdir) << " " << arma::norm(triads.slice(bp).col(2)) << std::endl;
            std::cout << initdir;
            std::cout << triads.slice(bp).col(2);
        }

		if (bp < num_bps) {
            T0 = *BPS[bp]->get_T0();
            Rz = Rotz(T0(2)+exc_tw_pbp);
            RzN = Rz*RzN;
        }
	}

	std::cout << "Initializing energies... ";
	init_energies();
    std::cout << "done \n";
    conf_initialized=true;
    return true;
}

bool Chain::gen_circular(unsigned number_bp, double supercoiling_density, const std::string& sequence) {
    std::cout << "Initializing circular chain configuration" << std::endl;
    if (conf_initialized) {
        std::cout << "Warning: Trying to generate circular configuration, but a configuration is already initialized" << std::endl;
        return false;
    }
    num_bp   = number_bp;
    num_bps  = number_bp;
    disc_len = data.get_disc_len();
    contour_len = num_bps*disc_len;
    set_seq(sequence,num_bp);
    set_closed_topology(true);
    config_type = "circular";

    // Initialize BPStep
    init_BPS();

    // Generate Configuration
    bp_pos = arma::zeros(3,num_bp);
    triads = arma::zeros(3,3,num_bp);

	double    theta = (2*M_PI)/num_bp;
	arma::mat Rot   = Rotz(theta);

	arma::mat triad = {{0.,1.,0.},{0.,0.,1.},{1.,0.,0.}};
	//~ triad = {{0.,0.,1.},{1.,0.,0.},{0.,1.,0.}};

	triads.slice(0) = triad;
	bp_pos.col(0)	= arma::zeros(3);

	for (unsigned bp=1;bp<num_bp;bp++) {
		bp_pos.col(bp) = bp_pos.col(bp-1)+disc_len*triad.col(2);
		triad = Rot*triad;
		triads.slice(bp) = triad;
	}

	/*
        Include intrinsic and excess twist
	*/
	dLK = sigma2dLk(supercoiling_density);
	double excess_twist_per_bps = (dLK*2*M_PI)/num_bps;

    double phi = 0;
	for (unsigned i=0;i<num_bp;i++) {
		Rot = Rotz(phi);
		triads.slice(i) = triads.slice(i)*Rot;
		phi = phi + (*BPS[i]->get_T0())(2) + excess_twist_per_bps;
	}

	init_energies();

    double check_dLK = 0;
    arma::colvec Theta;
    for (unsigned i=0;i<num_bps;i++) {
        Theta = *BPS[i]->get_Theta();
        check_dLK += Theta(2)/(2*M_PI);
    }
    std::cout << "check_dLK = " << check_dLK << std::endl;

    conf_initialized=true;
    return true;
}

bool Chain::gen_conf_from_restart(std::string restart_fn,int snapshot,std::string required_config_type) {

    std::cout << "Initializing chain from restart file" << std::endl;
    std::cout << " restart filename: " << restart_fn << std::endl;
    std::cout << " snapshot:         " << snapshot << std::endl;

    std::vector<RestartData> restarts = loadrestart(restart_fn);
    std::cout << "number of restarts: " << restarts.size() << std::endl;

    if (restarts.size()==0) {
        std::cout << "Error: Loading from restart file was unsuccessful!" << std::endl;
        std::exit(0);
        return false;
    }
    int selection=-1;
    if (snapshot != -1) {
        for (int i=0;i<restarts.size();i++) {
            if (restarts[i].snapshot == snapshot) {
                selection = i;
                break;
            }
        }
    }
    if (selection==-1) {
        selection = restarts.size() - 1;
    }

    if (restarts[selection].type == "circular") {

        num_bp      = restarts[selection].num_bp;
        num_bps     = num_bp;
        disc_len    = data.get_disc_len();
        contour_len = num_bps*disc_len;
        config_type = "circular";
        dLK         = restarts[selection].dLK;

        set_seq(restarts[selection].sequence,num_bp);
        set_closed_topology(true);

        bp_pos = restarts[selection].pos;
        triads = restarts[selection].triads;

        // Initialize BPStep
        init_BPS();
        init_energies();
        conf_initialized=true;
        return true;
    }
    if (restarts[selection].type == "linear") {

        num_bp      = restarts[selection].num_bp;
        num_bps     = num_bp-1;
        disc_len    = data.get_disc_len();
        contour_len = num_bps*disc_len;
        config_type = "linear";
        dLK         = restarts[selection].dLK;

        bp_pos = restarts[selection].pos;
        triads = restarts[selection].triads;

        set_seq(restarts[selection].sequence,num_bp);
        set_closed_topology(false);

        // Initialize BPStep
        init_BPS();
        init_energies();
        conf_initialized=true;
        return true;
    }

    std::cout << "Error: Loading from restart file was unsuccessful!" << std::endl;
    std::exit(0);

    return false;
}

bool Chain::init_custom_open_conf(arma::mat pos, arma::cube triads, const std::string& sequence)
{
    std::cout << "Initializing custom open chain configuration" << std::endl;
    if (conf_initialized) {
        std::cout << "Warning: Trying to initialize a custom open configuration, but a configuration is already initialized" << std::endl;
        return false;
    }

    this->bp_pos = pos;
    this->triads = triads;

    num_bp      = bp_pos.n_cols;
    num_bps     = num_bp-1;
    disc_len    = data.get_disc_len();
    /*
        check whether the given configuration is consistent with the specified discretization length
    */
    for (unsigned i=0;i<num_bps;i++) {
        if (std::abs(arma::norm(bp_pos.col(i+1)-bp_pos.col(i)) - disc_len ) > 1e-9) {
            std::cout << "Error: Chain::init_custom_open_conf: Given configuration is inconsistent with specified discretization length." << std::endl;
            std::cout << std::abs(arma::norm(bp_pos.col(i+1)-bp_pos.col(i)) - disc_len ) << std::endl;
            std::cout << arma::norm(bp_pos.col(i+1)-bp_pos.col(i)) << std::endl;
            std::cout << disc_len << std::endl;
            std::cout << i << std::endl;
            std::exit(0);
        }
    }
    std::cout << "seg: " << sequence << std::endl;

    contour_len = num_bps*disc_len;
    set_seq(sequence,num_bp);
    set_closed_topology(false);
    config_type = "custom_open";

    // Initialize BPStep
    init_BPS();
	std::cout << "Initializing energies... ";
	init_energies();
    std::cout << "done \n";

    dLK = cal_langowski_writhe_1a() + cal_twist();

    conf_initialized=true;
    return true;
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/*
    Set the orientation of the first and the last triad
*/

bool Chain::set_orientation_first_triad(arma::colvec& tangent) {
    arma::colvec tan, oldtan, rotax;
    arma::mat    rot_mat;
    double phi;
    arma::colvec oldpos,disp;

    tan = tangent/arma::norm(tangent);
    oldtan = triads.slice(0).col(2);

    if (std::abs(1-arma::dot(oldtan,tan)) < 1e-8) {
        return false;
    }
    if (std::abs(1+arma::dot(oldtan,tan)) < 1e-6) {
        rotax = triads.slice(0).col(1)*M_PI;
        rot_mat = getRotMat(rotax);
    }
    else {
        rotax = arma::cross(oldtan,tan);
        phi   = asin(arma::norm(rotax));
        std::cout << phi << std::endl;
        if (arma::dot(oldtan,tan)<0) {
            phi = M_PI-phi;
        }
        rotax    = phi*rotax/arma::norm(rotax);
        rot_mat = getRotMat(rotax);
    }
    triads.slice(0) = rot_mat*triads.slice(0);
    oldpos          = bp_pos.col(1);
    bp_pos.col(1)   = bp_pos.col(0) + triads.slice(0).col(2)*disc_len;
    disp            = bp_pos.col(1) - oldpos;

    for (unsigned i=2;i<num_bp;i++) {
        bp_pos.col(i) = bp_pos.col(i) + disp;
    }
    recal_energy();
    return true;
}

bool Chain::set_orientation_last_triad(arma::colvec& tangent) {
    arma::colvec tan, oldtan, rotax;
    arma::mat    rot_mat;
    double phi;

    tan = tangent/arma::norm(tangent);
    oldtan = triads.slice(num_bp-1).col(2);

    if (std::abs(1-arma::dot(oldtan,tan)) < 1e-8) {
        return false;
    }
    if (std::abs(1+arma::dot(oldtan,tan)) < 1e-6) {
        rotax = triads.slice(0).col(1)*M_PI;
        rot_mat = getRotMat(rotax);
    }
    else {
        rotax = arma::cross(oldtan,tan);
        phi   = asin(arma::norm(rotax));
        if (arma::dot(oldtan,tan)<0) {
            phi = M_PI-phi;
        }
        rotax   = phi*rotax/arma::norm(rotax);
        rot_mat = getRotMat(rotax);
    }

    triads.slice(num_bp-1) = rot_mat*triads.slice(num_bp-1);
    recal_energy();
    return true;
}



////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////// SET CONFIGURATION ///////////////////////////////////////////////////////////////////////////////

bool Chain::set_config(arma::mat* bp_pos, arma::cube* triads, bool closed, bool init) {
/*
    ToDo:
        dLK is not calculated in this function.
*/

    if (arma::size(*bp_pos)!=arma::size(this->bp_pos) || arma::size(*triads)!=arma::size(this->triads) ) {
        std::cout << "The new configuration has not the same size as the old one!" << std::endl;
        std::exit(0);
    }
    this->bp_pos  = *bp_pos;
    this->triads  = *triads;
    this->set_closed_topology(closed);
    if (init) {
        this->init_energies();
    }
    else {
        this->recal_energy();
    }
    return true;
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////// CALCULATE LINKING NUMBER  ///////////////////////////////////////////////////////////////////////////

double Chain::sigma2dLk(double sigma) {
    /*
        Calulates the excess linking number for the given supercoiling density sigma using the expression
        LK = LK_0(1+sigma)

        If the base pair steps include intrinsic twist LK_0 is deduced by the sum of the individual
        intrinsic twist contributions.

        If the steps do not include intrinsic twist, LK_0 is calculated based on the specified average
        intrinsic twist density (avg_intrinsic_twist_density), which is either preset or has been specified.
    */


    /*
        calculate LK_0 based on the intrinsic twist of the base pair steps
    */
    double LK_0=0;
    for (unsigned id=0;id<num_bps;id++) {
        LK_0+=(*BPS[id]->get_T0())(2);
    }
    LK_0/=2*M_PI;

    /*
        If this is (close to) zero, LK_0 is calculated based on the specified
        average intrisc twist density.
    */
    if (std::abs(LK_0)<1e-3) {
//        LK_0 = avg_intrinsic_twist_density*num_bps*disc_len/(2*M_PI);
        LK_0 = num_bps*disc_len/helical_repeat_length;
    }

    if (closed_topology) {
        double LK_0_mismatch=fpmod(LK_0,1);
        if (LK_0_mismatch > 0.5)   LK_0_mismatch=1-LK_0_mismatch;
        else                       LK_0_mismatch=-LK_0_mismatch;
        LK_0 = std::round(LK_0+LK_0_mismatch);
    }

    double LK = LK_0*(1+sigma);
    if (closed_topology) {
        double LK_mismatch=fpmod(LK,1);
        if (LK_mismatch > 0.5)   LK_mismatch=1-LK_mismatch;
        else                     LK_mismatch=-LK_mismatch;
        LK = std::round(LK+LK_mismatch);
    }

    return LK-LK_0;
}


////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///// CONFIGURATION INITIALIZATION AUXILIARY METHODS ///////////////////////////////////////////////////////////////

void Chain::init_BPS() {
    std::cout << "Creating BPStep objects... " << std::endl;
    // Init BPStep instances

    int next = 0;
    std::cout << "           |0%____________________________________________100%|" << std::endl;
    std::cout << " Progress:  " << std::flush;

    for (unsigned id=0;id<num_bps;id++) {
        BPStep * BPS_new = new BPStep(id, seq, data, disc_len,T, closed_topology);
        BPS_new->set_T0_subtract(T0_subtract);
        BPS.push_back(BPS_new);

        if (id>= 1.0*next/100.*num_bp) {
            std::cout << "\u2588" << std::flush;
            next += 2;
        }

    }
    std::cout << std::endl << "... done!" << std::endl;
    std::cout << num_bps << " BPS initialized\n";
    std::cout << "Initializing coupling neighbors... " << std::flush;
    // Init Neighbors
    for (unsigned id=0;id<num_bps;id++) {
        BPS[id]->init_neighbors(BPS);
    }
    std::cout << "done\n";
    std::cout << "Wire energies... " << std::flush;
    // Init Energies
    for (unsigned id=0;id<num_bps;id++) {
        BPS[id]->wire_energies();
    }
    std::cout << "done\n";
    std::cout << "Checking coupling consistency... \n";
    // Coupling Matrix Consistency Check
    bool consistent=true;
    for (unsigned id=0;id<num_bps;id++) {
        if (!BPS[id]->check_coupling_consistency((data.get_avg_inconsist() && AVG_COUP_MISSMATCH),AVG_COUP_HARMONIC)) {
            consistent = false;
        }
    }
    if (consistent) {
//        cout << ".--------------------------------------------." << std::endl;
//        cout << "| -> Coupling matrices are fully consistent. |" << std::endl;
//        cout << "'--------------------------------------------'" << std::endl;
        std::cout << " -> Coupling matrices are fully consistent." << std::endl;
    }
    else {
//        cout << ".------------------------------------------." << std::endl;
//        cout << "| -> Energetic Couplings are inconsistent! |" << std::endl;
//        cout << "'------------------------------------------'" << std::endl;
        std::cout << "  -> Energetic Couplings are inconsistent! " << std::endl;
        std::exit(0);
    }
    // calculate average stiffness and covariance matrix
    cal_avg_stiff(true);
}

void Chain::init_energies() {

    if (closed_topology) {
        for (unsigned bps=0;bps<num_bps-1;bps++) {
            BPS[bps]->propose_move(triads.slice(bps),triads.slice(bps+1));
        }
        BPS[num_bps-1]->propose_move(triads.slice(num_bps-1),triads.slice(0));
    }
    else {
        for (unsigned bps=0;bps<num_bps;bps++) {
            BPS[bps]->propose_move(triads.slice(bps),triads.slice(bps+1));
        }
    }
    for (unsigned bps=0;bps<num_bps;bps++) {
        BPS[bps]->eval_delta_energy();
    }
    for (unsigned bps=0;bps<num_bps;bps++) {
        BPS[bps]->set_move(true);
    }
}

void Chain::set_seq(const std::string& sequence, int num_bp) {
/*
Sets the sequence of the chain. If num_bp is a multiple of the length of
the sequence a repetitive sequence is assigned.
*/
    unsigned lenseq = sequence.length();
    if ((int)lenseq==num_bp) { seq = sequence;}
    else {
        if (num_bp%lenseq==0) {
            int repeats=num_bp/lenseq;
            seq = "";
            for (int i=0;i<repeats;i++) { seq = seq+sequence;}
        }
        else {
            throw std::domain_error("Chain::set_seq(): Provided sequence cannot be matched with set length.");
        }
    }
}


////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////// CONFIGURATION MODIFICATIONS  /////////////////////////////////////////////////////////////////////////

void Chain::set_Delta_Lk(double new_dLK) {
    if (closed_topology) {
        new_dLK = std::round(new_dLK);
    }
    double change_dLk = new_dLK-dLK;
    double change_twist_per_bps = (change_dLk*2*M_PI)/(num_bps);
    double phi = change_twist_per_bps;
    arma::mat Rot;
	for (unsigned i=1;i<num_bp;i++) {
		Rot = Rotz(phi);
		triads.slice(i) = triads.slice(i)*Rot;
		phi = phi + change_twist_per_bps;
	}
    init_energies();
    dLK = new_dLK;

    if (!check_link_conservation()) {
        std::cout << "Error occured during change of linking number!" << std::endl;
    }
}


void Chain::set_Delta_Lk(double new_dLK, int id_first, int id_last) {
    if (closed_topology) {
        if (id_first == 0 && id_last == num_bp-1) {
            set_Delta_Lk(new_dLK);
        }
        double ref_new_dLK = new_dLK;
        new_dLK = int(new_dLK);
        if (std::abs(ref_new_dLK - new_dLK) > 1e-8) {
            std::cout << "Warning: The linking number of a closed chain has to be an integer. The provided value was typecast to integer." << std::endl;
        }
        throw std::domain_error("Chain::set_Delta_Lk(): Range selection for setting Linking number not yet supported for closed chains.");
        return;
    }
    if (id_last < id_first || id_first < 0 || id_last >= num_bp) {
        throw std::domain_error("Chain::set_Delta_Lk(): Invalid range for setting linking number. ");
        id_last = num_bp-1;
    }

    int num_changed_bps         = id_last-id_first;
    double change_dLk           = new_dLK-dLK;
    double change_twist_per_bps = (change_dLk*2*M_PI)/(num_changed_bps);

    double phi = change_twist_per_bps;
    arma::mat Rot;
	for (unsigned i=id_first+1;i<=id_last;i++) {
		Rot = Rotz(phi);
		triads.slice(i) = triads.slice(i)*Rot;
		phi = phi + change_twist_per_bps;
	}
	for (unsigned i=id_last+1;i<num_bp;i++) {
        triads.slice(i) = triads.slice(i)*Rot;
	}

    init_energies();
    dLK = new_dLK;

    if (!check_link_conservation()) {
        std::cout << "Error occured during change of linking number!" << std::endl;
    }
}


void Chain::set_sigma(double sigma) {
    set_Delta_Lk(sigma2dLk(sigma));
}

void Chain::set_sigma(double sigma, int id_first, int id_last) {
    set_Delta_Lk(sigma2dLk(sigma),id_first,id_last);
}
