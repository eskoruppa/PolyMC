#include "MCS_Kinkxy.h"

/*
TODO:
create list of kinkxy bps in constructor and only select hinges from those elements
*/


MCS_Kinkxy::MCS_Kinkxy(Chain * ch,const std::vector<long long int> & seedseq)
: MCStep(ch,seedseq)
{
    move_name="MCS_Kinkxy";
    // initialize range of random hinge selector
    decltype(genhinge.param()) new_range(1, num_bp-1);
    genhinge.param(new_range);

//    sigfac = 4.130;
//    fac_sigma = 4;
//
    arma::mat covmat = chain->get_avg_cov();
//    sigma = fac_sigma*1./3*(sqrt(covmat(0,0))+sqrt(covmat(1,1))+sqrt(covmat(2,2)));

    double lb = 2.*disc_len/(covmat(0,0) + covmat(1,1));
    double f = chain->get_force();
    double kT    = chain->get_kT();
    double tail_length;

//    double lt = disc_len/covmat(2,2);
//    std::cout << lb << " " << lt << std::endl;
//    std::exit(0);

    if (f==0) {
        tail_length = lb;
    }
    else {
        double xi = std::sqrt(kT*lb/f);
        if (xi < lb/10) {
            tail_length = lb/10;
        }
        else if (xi > lb) {
            tail_length = lb;
        }
        else {
            tail_length = xi;
        }
    }
    num_in_tail = int(std::ceil(tail_length/disc_len));

    if (num_bp < 4*num_in_tail) {
        std::cout << "\nError: MCStep MCS_Kinkxy requires a chain length of at least 4*lb" << std::endl;
        std::cout << "  bp required:  " << num_in_tail*4 << std::endl;
        std::cout << "  bp contained: " << num_bp << std::endl;
        std::exit(0);
    }

//    std::cout << lb << std::endl;
//    std::cout << disc_len/covmat(2,2) << std::endl;
//    std::cout << num_in_tail << std::endl;
//    std::exit(0);


    // intervals of moved segments for exluded volume interactions
    moved_intervals.push_back({0,-1,0});
    moved_intervals.push_back({-1,-1,2000});
    moved_intervals.push_back({-1,num_bp-1,1});

    suitablility_key[MCS_FORCE_ACTIVE]              = true;  // force_active
    suitablility_key[MCS_TORQUE_ACTIVE]             = true;  // torque_active
    suitablility_key[MCS_TOPOLOGY_CLOSED]           = false; // topology_closed
    suitablility_key[MCS_FIXED_LINK]                = true;  // link_fixed
    suitablility_key[MCS_FIXED_TERMINI]             = true;  // fixed_termini
    suitablility_key[MCS_FIXED_TERMINI_RADIAL]      = true;  // fixed_termini_radial
    suitablility_key[MCS_FIXED_FIRST_ORIENTATION]   = true;  // fixed_first_orientation
    suitablility_key[MCS_FIXED_LAST_ORIENTATION]    = true;  // fixed_last_orientation

    requires_EV_check=true;
}


MCS_Kinkxy::~MCS_Kinkxy() {

}

void MCS_Kinkxy::update_settings() {
    arma::mat covmat = chain->get_avg_cov();
    sigma = fac_sigma*1./3*(sqrt(covmat(0,0))+sqrt(covmat(1,1))+sqrt(covmat(2,2)));
}

bool MCS_Kinkxy::MC_move() {
    int    hingeID  = genhinge(gen);
    changed_bps = {};

    // if the hinge BPS is not of type kinkxy this move does nothing
    if (*BPS[hingeID-1]->get_method_diag() != "kinkxy") {
        // will become redundant once hinges are only chosen from valid list
        return false;
    }

    if (hingeID <= num_in_tail) {
        return right_tail_move(hingeID);
    }
    if (num_bp-hingeID-1 <= num_in_tail) {
        return left_tail_move(hingeID);
    }

    return double_tail_move(hingeID);
}


///////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////
// Double Tail Move

bool MCS_Kinkxy::double_tail_move(int idB) {

    int idBm1 = idB-1;
    int idA   = idB-num_in_tail;
    int idC   = idB+num_in_tail;

    int idAm1 = idA-1;
    int idBp1 = idB+1;
    int idCp1 = idC+1;

    double choice   = uniformdist(gen);
    double deltaE   = 0;

    arma::mat Ta,Tb;
    arma::colvec current_Theta, swapped_Theta;
    arma::mat Rkink,Rsq;

    arma::colvec    phi,phi_step;
    arma::mat       R,Rstep,RstepT;
    arma::mat       Raccu;

    arma::colvec dRi,dRf;

    // triads at kink
    Ta = triads->slice(idBm1);
    Tb = triads->slice(idB);

    // get current Theta at kink
    current_Theta = BPS[idBm1]->Triads2Theta(Ta,Tb);

    // get kink swapped Theta if such can be found
    if (!get_kink_swap_Theta(idBm1, current_Theta, swapped_Theta)) {
        // no valid move found (return true without changing configuration to avoid reverting to backup config)
        return false;
    }

    // initial Vector between R_idA and R_idC
    dRi = pos->col(idC) - pos->col(idA);

    // Rotation matrix relating transformed Ta and Tb in frame of Ta
    Rkink = BPS[idBm1]->Theta2Rotmat(swapped_Theta);

    // R^2
    Rsq = Ta*Rkink*Tb.t();

    // rotation vectors for half rotation, full rotation and tail step
    phi         = ExtractTheta(Rsq)/2;
    phi_step    = phi/num_in_tail;

    // rotation matrices for half rotation, full rotation and tail step
    R           = getRotMat(phi);
    Rstep       = getRotMat(phi_step);
    RstepT      = Rstep.t();
//    RstepT      = getRotMat(-phi_step);

    // create configuration backup
    set_trial_backup();

    // rotate left tail
    Raccu = arma::eye(3,3);
    for (int i=idA;i<idB;i++) {
        Raccu *= RstepT;
        triads->slice(i)   = Raccu * triads->slice(i);
        pos->col(i+1)      = pos->col(i) + triads->slice(i).col(2)*disc_len;
        BPS[i-1]->propose_move(triads->slice(i-1),triads->slice(i));
    }

    // rotate kink and right tail
    Raccu = R;
    for (int i=idB;i<idC;i++) {
        triads->slice(i) = Raccu * triads->slice(i);
        pos->col(i+1)      = pos->col(i) + triads->slice(i).col(2)*disc_len;
        BPS[i-1]->propose_move(triads->slice(i-1),triads->slice(i));
        Raccu *= RstepT;
    }
    BPS[idC-1]->propose_move(triads->slice(idC-1),triads->slice(idC));

    // final Vector between R_idA and R_idC
    dRf = pos->col(idC) - pos->col(idA);

    // Evaluate elastic energy
    for (int i=idAm1;i<idC;i++) {
        deltaE += BPS[i]->eval_delta_energy();
    }
    changed_bps = arma::zeros<arma::ivec>(idC-idAm1);
    for (int i=0;i<idC-idAm1;i++) {
        changed_bps(i) = i+idAm1;
    }

    // Evaluate change in work
	if (chain->force_active() == true) {
        deltaE += arma::dot(chain->get_beta_force_vec(),dRi-dRf);
	}

    // Metropolis step
	if (exp(-deltaE) <= choice) {
        revert_to_trial_backup();
		return false;
    }
	else {
        // translate segments to the right of idC
        for (int i=idC+1;i<num_bp;i++) {
            pos->col(i)      = pos->col(i-1) + triads->slice(i-1).col(2)*disc_len;
        }

        // the first one is permanently set to 0
        moved_intervals[0](EV_TO)   = idA;
        moved_intervals[1](EV_FROM) = idA+1;
        moved_intervals[1](EV_TO)   = idC;
        moved_intervals[2](EV_FROM) = idC+1;
        // the last one is permanently set to num_bp-1
    }

    return true;
}

///////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////
// Left Tail Move

bool MCS_Kinkxy::left_tail_move(int idB) {

    int idA   = idB - 2*num_in_tail;
    int idBm1 = idB - 1;
    int idAm1 = idA - 1;

    double choice   = uniformdist(gen);
    double deltaE   = 0;

    arma::mat Ta,Tb;
    arma::colvec current_Theta, swapped_Theta;
    arma::mat Rkink;

    arma::colvec    phi,phi_step;
    arma::mat       R,Rstep;
    arma::mat       Raccu;

    arma::colvec dRi,dRf;

    // triads at kink
    Ta = triads->slice(idBm1);
    Tb = triads->slice(idB);

    // get current Theta at kink
    current_Theta = BPS[idBm1]->Triads2Theta(Ta,Tb);

    // get kink swapped Theta if such can be found
    if (!get_kink_swap_Theta(idBm1, current_Theta, swapped_Theta)) {
        // no valid move found (return true without changing configuration to avoid reverting to backup config)
        return false;
    }

    // initial Vector between R_idA and R_idC
    dRi = pos->col(idB) - pos->col(idA);

    // Rotation matrix relating transformed Ta and Tb in frame of Ta
    Rkink = BPS[idBm1]->Theta2Rotmat(swapped_Theta);

    // Lab Frame rotation matrix inducing kinkswap between Ta and Tb
    R = Tb*Rkink.t()*Ta.t();

    // rotation vectors for full rotation and tail step
    phi         = ExtractTheta(R);
    phi_step    = phi/(num_in_tail*2);

    // rotation matrices for tail step
    R           = getRotMat(phi);
    Rstep       = getRotMat(phi_step);

    // create configuration backup
    set_trial_backup();

    Raccu = arma::eye(3,3);
    for (int i=idA;i<idB;i++) {
        Raccu *= Rstep;
        triads->slice(i) = Raccu * triads->slice(i);
        pos->col(i+1)      = pos->col(i) + triads->slice(i).col(2)*disc_len;
        BPS[i-1]->propose_move(triads->slice(i-1),triads->slice(i));
    }
    BPS[idBm1]->propose_move(triads->slice(idBm1),triads->slice(idB));

    // final Vector between R_idA and R_idC
    dRf = pos->col(idB) - pos->col(idA);

    // Evaluate elastic energy
    for (int i=idAm1;i<idB;i++) {
        deltaE += BPS[i]->eval_delta_energy();
    }
    changed_bps = arma::zeros<arma::ivec>(idB-idAm1);
    for (int i=0;i<idB-idAm1;i++) {
        changed_bps(i) = i+idAm1;
    }

    // Evaluate change in work
	if (chain->force_active() == true) {
        deltaE += arma::dot(chain->get_beta_force_vec(),dRi-dRf);
//        std::cout << arma::dot(chain->get_beta_force_vec(),dRi-dRf) << std::endl;
	}

    // Metropolis step
	if (exp(-deltaE) <= choice) {
        revert_to_trial_backup();
		return false;
    }
	else {
        // translate segments to the right of idC
        for (int i=idB+1;i<num_bp;i++) {
            pos->col(i)      = pos->col(i-1) + triads->slice(i-1).col(2)*disc_len;
        }

        // the first one is permanently set to 0
        moved_intervals[0](EV_TO)   = idA;
        moved_intervals[1](EV_FROM) = idA+1;
        moved_intervals[1](EV_TO)   = idB;
        moved_intervals[2](EV_FROM) = idB+1;
        // the last one is permanently set to num_bp-1
    }
    return true;
}

///////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////
// Right Tail Move

bool MCS_Kinkxy::right_tail_move(int idA) {

    int idAm1 = idA - 1;
    int idB   = idA + 2*num_in_tail;

    double choice   = uniformdist(gen);
    double deltaE   = 0;

    arma::mat Ta,Tb;
    arma::colvec current_Theta, swapped_Theta;
    arma::mat Rkink;

    arma::colvec    phi,phi_step;
    arma::mat       R,Rstep,RstepT;
    arma::mat       Raccu;

    arma::colvec dRi,dRf;

    // triads at kink
    Ta = triads->slice(idAm1);
    Tb = triads->slice(idA);

    // get current Theta at kink
    current_Theta = BPS[idAm1]->Triads2Theta(Ta,Tb);

    // get kink swapped Theta if such can be found
    if (!get_kink_swap_Theta(idAm1, current_Theta, swapped_Theta)) {
        // no valid move found (return true without changing configuration to avoid reverting to backup config)
        return false;
    }

    // initial Vector between R_idA and R_idC
    dRi = pos->col(idB) - pos->col(idA);

    // Rotation matrix relating transformed Ta and Tb in frame of Ta
    Rkink = BPS[idAm1]->Theta2Rotmat(swapped_Theta);

    // Lab Frame rotation matrix inducing kinkswap between Ta and Tb
    R = Ta*Rkink*Tb.t();

    // rotation vectors for full rotation and tail step
    phi         = ExtractTheta(R);
    phi_step    = phi/(num_in_tail*2);

    // rotation matrices for tail step
    R           = getRotMat(phi);
    Rstep       = getRotMat(phi_step);
    RstepT      = Rstep.t();

    // create configuration backup
    set_trial_backup();

    // rotate kink
    triads->slice(idA) = R * triads->slice(idA);
    pos->col(idA+1)    = pos->col(idA) + triads->slice(idA).col(2)*disc_len;
    BPS[idA-1]->propose_move(triads->slice(idA-1),triads->slice(idA));

    // rotate kink and right tail
    Raccu = R;
    for (int i=idA+1;i<idB;i++) {
        Raccu *= RstepT;
        triads->slice(i) = Raccu * triads->slice(i);
        pos->col(i+1)      = pos->col(i) + triads->slice(i).col(2)*disc_len;
        BPS[i-1]->propose_move(triads->slice(i-1),triads->slice(i));
    }
    BPS[idB-1]->propose_move(triads->slice(idB-1),triads->slice(idB));

    // final Vector between R_idA and R_idC
    dRf = pos->col(idB) - pos->col(idA);

    // Evaluate elastic energy
    for (int i=idAm1;i<idB;i++) {
        deltaE += BPS[i]->eval_delta_energy();
    }
    changed_bps = arma::zeros<arma::ivec>(idB-idAm1);
    for (int i=0;i<idB-idAm1;i++) {
        changed_bps(i) = i+idAm1;
    }

    // Evaluate change in work
	if (chain->force_active() == true) {
        deltaE += arma::dot(chain->get_beta_force_vec(),dRi-dRf);
//        std::cout << arma::dot(chain->get_beta_force_vec(),dRi-dRf) << std::endl;
	}

    // Metropolis step
	if (exp(-deltaE) <= choice) {
        revert_to_trial_backup();
		return false;
    }
	else {
        // translate segments to the right of idC
        for (int i=idB+1;i<num_bp;i++) {
            pos->col(i)      = pos->col(i-1) + triads->slice(i-1).col(2)*disc_len;
        }

        // the first one is permanently set to 0
        moved_intervals[0](EV_TO)   = idA;
        moved_intervals[1](EV_FROM) = idA+1;
        moved_intervals[1](EV_TO)   = idB;
        moved_intervals[2](EV_FROM) = idB+1;
        // the last one is permanently set to num_bp-1
    }
    return true;
}



///////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////
// Additional methods


bool MCS_Kinkxy::get_kink_swap_Theta(int bpsid, const arma::colvec & Theta, arma::colvec & newTheta) {
    std::vector<double> attr = *BPS[bpsid]->get_EvalEnergy_diag()->get_attributes();

    if (attr[KINKXY_left_active] == 0 && attr[KINKXY_right_active] == 0) {
        return false;
    }

    int kink_state = get_kink_state(Theta ,attr);

    // in left kinked state -> swap to neutral
    if (kink_state == -1) {
        return get_theta_swap_left_kink(Theta,attr,newTheta);
    }

    // in right kinked state -> swap to neutral
    if (kink_state == 1) {
        return get_theta_swap_right_kink(Theta,attr,newTheta);
    }

    // in neutral state
    if (attr[KINKXY_left_active] == 1) {
        if (attr[KINKXY_right_active] == 1) {
            // if both active chose on kink at random
            if (randbinary(gen) == 0) {
                return get_theta_swap_left_kink(Theta,attr,newTheta);
            }
            else {
                return get_theta_swap_right_kink(Theta,attr,newTheta);
            }
        }
        return get_theta_swap_left_kink(Theta,attr,newTheta);
    }
    else {
        return get_theta_swap_right_kink(Theta,attr,newTheta);
    }
}

int MCS_Kinkxy::get_kink_state(const arma::colvec & Theta ,const std::vector<double> & attr) {
/*
    returns length 1 vector indicating whether bps is in kinked state
    unkinked:  0
    left:     -1
    right:    +1
*/
    double x;
    x =  attr[KINKXY_cosTk] * Theta(1) - attr[KINKXY_sinTk] * Theta(2);
    // x contribution
    if (attr[KINKXY_left_active] == 1 && x<attr[KINKXY_xcl]) {
        return -1;
    }
    if (attr[KINKXY_right_active] == 1 && x>attr[KINKXY_xcr]) {
        return 1;
    }
    return 0;
}

bool MCS_Kinkxy::get_theta_swap_left_kink(const arma::colvec & Theta, const std::vector<double> & attr, arma::colvec & newTheta) {
    if (attr[KINKXY_left_active] == 0) {
        // if left kink is not active return current Theta
        return false;
    }

    double x,y;
    x =  attr[KINKXY_cosTk]* Theta(1) - attr[KINKXY_sinTk]* Theta(2);
    y =  attr[KINKXY_sinTk]* Theta(1) + attr[KINKXY_cosTk]* Theta(2);

    if (attr[KINKXY_right_active] == 1 && x>attr[KINKXY_xcr]) {
        // if BPS is currently in right kink return current Theta (no move)
        return false;
    }

    double x_new;
    if (x<attr[KINKXY_xcl]) {
        // is currently in left kink
        x_new = x - attr[KINKXY_xl];
    }
    else {
        // is currently not kinked
        x_new = x + attr[KINKXY_xl];
    }

    newTheta = arma::zeros(3);
    newTheta(0) = Theta(0);
    newTheta(1) =  attr[KINKXY_cosTk] * x_new + attr[KINKXY_sinTk] * y;
    newTheta(2) = -attr[KINKXY_sinTk] * x_new + attr[KINKXY_cosTk] * y;
    return true;
}

bool MCS_Kinkxy::get_theta_swap_right_kink(const arma::colvec & Theta, const std::vector<double> & attr, arma::colvec & newTheta) {
    if (attr[KINKXY_right_active] == 0) {
        // if right kink is not active return current Theta
        return false;
    }

    double x,y;
    x =  attr[KINKXY_cosTk] * Theta(1) - attr[KINKXY_sinTk] * Theta(2);
    y =  attr[KINKXY_sinTk] * Theta(1) + attr[KINKXY_cosTk] * Theta(2);

    if (attr[KINKXY_left_active] == 1 && x<attr[KINKXY_xcl]) {
        // if BPS is currently in left kink return current Theta (no move)
        return false;
    }

    double x_new;
    if (x>attr[KINKXY_xcr]) {
        // is currently in right kink
        x_new = x - attr[KINKXY_xr];
    }
    else {
        // is currently not kinked
        x_new = x + attr[KINKXY_xr];
    }

    newTheta = arma::zeros(3);
    newTheta(0) = Theta(0);
    newTheta(1) =   attr[KINKXY_cosTk] * x_new + attr[KINKXY_sinTk] * y;
    newTheta(2) = - attr[KINKXY_sinTk] * x_new + attr[KINKXY_cosTk] * y;
    return true;
}


/////////////////////////////////////////////////
// Test reversibility of move
bool MCS_Kinkxy::test_reverse(int idB) {

    int idBm1 = idB-1;
    int idA   = idB-num_in_tail;
    int idC   = idB+num_in_tail;

    int idAm1 = idA-1;
    int idBp1 = idB+1;
    int idCp1 = idC+1;

    double choice   = uniformdist(gen);
    double deltaE   = 0;

    arma::mat Ta,Tb;
    arma::colvec current_Theta, swapped_Theta;
    arma::mat Rkink,Rsq;

    arma::colvec    phi,phi_step;
    arma::mat       R,Rstep,RstepT;
    arma::mat       Raccu;

    arma::colvec dRi,dRf;

    // triads at kink
    Ta = triads->slice(idBm1);
    Tb = triads->slice(idB);

    // get current Theta at kink
    current_Theta = BPS[idBm1]->Triads2Theta(Ta,Tb);

    // get kink swapped Theta if such can be found
    if (!get_kink_swap_Theta(idBm1, current_Theta, swapped_Theta)) {
        // no valid move found (return true without changing configuration to avoid reverting to backup config)
        return true;
    }

    // initial Vector between R_idA and R_idC
    dRi = pos->col(idC) - pos->col(idA);

    // Rotation matrix relating transformed Ta and Tb in frame of Ta
    Rkink = BPS[idBm1]->Theta2Rotmat(swapped_Theta);

    // R^2
    Rsq = Ta*Rkink*Tb.t();

    // rotation vectors for half rotation, full rotation and tail step
    phi         = ExtractTheta(Rsq)/2;
    phi_step    = phi/num_in_tail;

    // rotation matrices for half rotation, full rotation and tail step
    R           = getRotMat(phi);
    Rstep       = getRotMat(phi_step);
    RstepT      = Rstep.t();
//    RstepT      = getRotMat(-phi_step);

    // create configuration backup
    set_trial_backup();

    // rotate left tail
    Raccu = arma::eye(3,3);
    for (int i=idA;i<idB;i++) {
        Raccu *= RstepT;
        triads->slice(i) = Raccu * triads->slice(i);
        pos->col(i+1)      = pos->col(i) + triads->slice(i).col(2)*disc_len;
        BPS[i-1]->propose_move(triads->slice(i-1),triads->slice(i));
    }

    // rotate kink and right tail
    Raccu = R;
    for (int i=idB;i<idC;i++) {
        triads->slice(i) = Raccu * triads->slice(i);
        pos->col(i+1)      = pos->col(i) + triads->slice(i).col(2)*disc_len;
        BPS[i-1]->propose_move(triads->slice(i-1),triads->slice(i));
        Raccu *= RstepT;
    }
    BPS[idC-1]->propose_move(triads->slice(idC-1),triads->slice(idC));

    // final Vector between R_idA and R_idC
    dRf = pos->col(idC) - pos->col(idA);

    // Evaluate elastic energy
    for (int i=idAm1;i<idC;i++) {
        deltaE += BPS[i]->eval_delta_energy();
    }
    changed_bps = arma::zeros<arma::ivec>(idC-idAm1);
    for (int i=0;i<idC-idAm1;i++) {
        changed_bps(i) = i+idAm1;
    }

    // Evaluate change in work
//	if (chain->force_active() == true) {
//        deltaE += arma::dot(chain->get_beta_force_vec(),dRi-dRf);
////        std::cout << arma::dot(chain->get_beta_force_vec(),dRi-dRf) << std::endl;
//	}

//	std::cout << "dE      = " << deltaE << std::endl;

    // Metropolis step
//	if (exp(-deltaE) <= choice) {
//        revert_to_trial_backup();
//		return false;
//    }
	if (true) {
        // translate segments to the right of idC
        for (int i=idC+1;i<num_bp;i++) {
            pos->col(i)      = pos->col(i-1) + triads->slice(i-1).col(2)*disc_len;
        }

        // the first one is permanently set to 0
        moved_intervals[0](EV_TO)   = idA;
        moved_intervals[1](EV_FROM) = idA+1;
        moved_intervals[1](EV_TO)   = idC;
        moved_intervals[2](EV_FROM) = idC+1;
        // the last one is permanently set to num_bp-1
    }

    return true;
}



