#include "MCS_TopolPerm.h"

/*
    Topology Permutation Move
*/


MCS_TopolPerm::MCS_TopolPerm(Chain * ch,const std::vector<long long int> & seedseq, int trigger_every)
: MCStep(ch,seedseq),
trigger_every(trigger_every)
{
    move_name = "MCS_TopolPerm";
    /*
    Initializes the controling variables of the MCS_TopolPerm move
    */

    closed = chain->topology_closed();
    /*
        Quick evaluation is posible for local isotropic hamiltonians. Requires elastic energy evaluation
        only at the hinge site and allows for triad rotations to be carried out only if move is accepted.
    */

    nbpm1 = num_bp - 1;
    decltype(direction_selection.param()) new_range(0, 1);
    direction_selection.param(new_range);

    // gen_full_trial_conf = true;
    changed_bps = arma::regspace<arma::ivec>(0, 1, num_bp-2);

    suitablility_key[MCS_FORCE_ACTIVE]              = true;  // force_active
    suitablility_key[MCS_TORQUE_ACTIVE]             = true;  // torque_active
    suitablility_key[MCS_TOPOLOGY_CLOSED]           = true;  // topology_closed
    suitablility_key[MCS_FIXED_LINK]                = false;  // link_fixed
    suitablility_key[MCS_FIXED_TERMINI]             = true;  // fixed_termini
    suitablility_key[MCS_FIXED_TERMINI_RADIAL]      = true;  // fixed_termini_radial
    suitablility_key[MCS_FIXED_FIRST_ORIENTATION]   = true;  // fixed_first_orientation
    suitablility_key[MCS_FIXED_LAST_ORIENTATION]    = true;  // fixed_last_orientation

    requires_EV_check=false;
}


MCS_TopolPerm::~MCS_TopolPerm() {

}

void MCS_TopolPerm::update_settings() {

}

bool MCS_TopolPerm::MC_move() {

    if (count_step % trigger_every != 0) {
        return false;
    }

    double change_turns = 0;
    if (direction_selection(gen) == 0) {
        change_turns = -1;
    }
    else {
        change_turns = 1;
    }

    double change_twist_per_bps = (change_turns*2*M_PI)/(nbpm1);
    double phi;
    arma::mat Rot;

    arma::cube triads_cp = *triads;
    for (unsigned i=0;i<nbpm1;i++) {
        // std::cout << i << "/" << num_bp << std::endl;
        phi = (i+1)*change_twist_per_bps;
        Rot = Rotz(phi);
        triads_cp.slice(i+1) = triads_cp.slice(i+1)*Rot;
        BPS[i]->propose_move(triads_cp.slice(i), triads_cp.slice(i+1));
    }
    double deltaE = 0;  
    for (unsigned i=0;i<nbpm1;i++) {
        deltaE += BPS[i]->eval_delta_energy();
    }
    
    if (exp(-deltaE) <= uniformdist(gen)) {
        return false;
    }
    
    // std::cout << change_turns << " " << count_step << std::endl;
    // std::cout << "deltaE = " << deltaE << std::endl;
    // std::cout << "Tw = " << chain->cal_twist(0,num_bp) << std::endl;
    *triads = triads_cp;

    // for (unsigned i=0;i<changed_bps.n_elem;i++) {
    //     BPS[changed_bps(i)]->set_move(true);
    // }
    // std::cout << "Tw = " << chain->cal_twist(0,num_bp) << std::endl;
    return true;

}
