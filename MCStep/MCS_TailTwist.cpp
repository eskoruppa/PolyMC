#include "MCS_TailTwist.h"


MCS_TailTwist::MCS_TailTwist(Chain * ch,const std::vector<long long int> & seedseq)
: MCStep(ch,seedseq)
{
    move_name = "MCS_TailTwist";

    sigfac = 2;
    if (chain->torsional_trap_active()) {
        sigma = sigfac*std::abs(chain->get_torstrap_dLK_aim())*2*M_PI;
    }
    else {
        arma::mat covmat = chain->get_avg_cov();
        sigma = sigfac*sqrt(covmat(2,2));
    }

    // initialize range of random hinge selector
    decltype(genhinge.param()) new_range(1, num_bp-1);
    genhinge.param(new_range);

    /*
        Quick evaluation is posible for local isotropic hamiltonians. Requires elastic energy evaluation
        only at the hinge site and allows for triad rotations to be carried out only if move is accepted.
    */
    if (chain->local_isotropic_hamiltonian()) {
        quick_eval = true;
    }
    else {
        quick_eval = false;
    }

    suitablility_key[MCS_FORCE_ACTIVE]              = true;  // force_active
    suitablility_key[MCS_TORQUE_ACTIVE]             = true;  // torque_active
    suitablility_key[MCS_TOPOLOGY_CLOSED]           = false;  // topology_closed
    suitablility_key[MCS_FIXED_LINK]                = false;  // link_fixed
    suitablility_key[MCS_FIXED_TERMINI]             = false;  // fixed_termini
    suitablility_key[MCS_FIXED_TERMINI_RADIAL]      = true;  // fixed_termini_radial
    suitablility_key[MCS_FIXED_FIRST_ORIENTATION]   = true;  // fixed_first_orientation
    suitablility_key[MCS_FIXED_LAST_ORIENTATION]    = false;  // fixed_last_orientation

    requires_EV_check=false;
}

MCS_TailTwist::MCS_TailTwist(Chain * ch,const std::vector<long long int> & seedseq, int range_id1)
: MCStep(ch,seedseq)
{
    move_name = "MCS_TailTwist";

    sigfac = 2;
    if (chain->torsional_trap_active()) {
        sigma = sigfac*std::abs(chain->get_torstrap_dLK_aim())*2*M_PI;
    }
    else {
        arma::mat covmat = chain->get_avg_cov();
        sigma = sigfac*sqrt(covmat(2,2));
    }

    // initialize range of random hinge selector
    if (range_id1 < 0 || range_id1 > num_bp-1) {
        std::cout << "Error invalid Range selection in MCS_TailTwist. range_id1 needs to be within the range [1,num_bp-1]!" << std::endl;
        std::exit(0);
    }
    decltype(genhinge.param()) new_range(range_id1, num_bp-1);
    genhinge.param(new_range);

    /*
        Quick evaluation is posible for local isotropic hamiltonians. Requires elastic energy evaluation
        only at the hinge site and allows for triad rotations to be carried out only if move is accepted.
    */
    if (chain->local_isotropic_hamiltonian()) {
        quick_eval = true;
    }
    else {
        quick_eval = false;
    }

    suitablility_key[MCS_FORCE_ACTIVE]              = true;  // force_active
    suitablility_key[MCS_TORQUE_ACTIVE]             = true;  // torque_active
    suitablility_key[MCS_TOPOLOGY_CLOSED]           = false;  // topology_closed
    suitablility_key[MCS_FIXED_LINK]                = false;  // link_fixed
    suitablility_key[MCS_FIXED_TERMINI]             = false;  // fixed_termini
    suitablility_key[MCS_FIXED_TERMINI_RADIAL]      = true;  // fixed_termini_radial
    suitablility_key[MCS_FIXED_FIRST_ORIENTATION]   = true;  // fixed_first_orientation
    suitablility_key[MCS_FIXED_LAST_ORIENTATION]    = false;  // fixed_last_orientation

    requires_EV_check=false;
}


MCS_TailTwist::~MCS_TailTwist() {

}

void MCS_TailTwist::update_settings() {
    if (chain->torsional_trap_active()) {
        sigma = sigfac*std::abs(chain->get_torstrap_dLK_aim())*2*M_PI;
    }
    else {
        arma::mat covmat = chain->get_avg_cov();
        sigma = sigfac*sqrt(covmat(2,2));
    }
}

bool MCS_TailTwist::MC_move() {
    int    hinge_id = genhinge(gen);
    int    idm1     = hinge_id-1;
    double choice   = uniformdist(gen);
    double deltaE,dtheta;

    dtheta        = normaldist(gen)*sigma;
    arma::mat Rz  = Rotz(dtheta);

    if (quick_eval) {
        /*
            If the hamiltonian is local and bending is isotropic the energy will only be evaluated at the hinge
            and at the terminal bead if the torsional trap or an external torque constraint is active.
        */

        arma::mat T_hinge, T_last;
        T_hinge = triads->slice(hinge_id)*Rz;
        BPS[idm1]->propose_move(triads->slice(idm1), T_hinge);
        deltaE = BPS[idm1]->eval_delta_energy();
        changed_bps = {idm1};

        /*
            Contribution due to torque or torsional trap
        */
        deltaE += chain->propose_terminus_rotation( triads->slice(num_bp-1) * Rz, dtheta );

//        std::cout << sigma << std::endl;
//        std::cout << chain->get_torstrap_dLK_aim() << std::endl;
//        std::cout << Rz << std::endl;
//        std::cout << deltaE << std::endl;

        if (exp(-deltaE) <= choice) {
            chain->set_terminus_rotation(false);
            return false;
        }

        triads->slice(hinge_id) = T_hinge;
        for (int i=hinge_id+1;i<num_bp;i++) {
            triads->slice(i) = triads->slice(i)*Rz;
        }
        chain->set_terminus_rotation(true);
        return true;
    }
    else {
        arma::cube triads_cp = *triads;
        for (int i=hinge_id;i<num_bp;i++) {
            triads_cp.slice(i) = triads_cp.slice(i)*Rz;
            BPS[i-1]->propose_move(triads_cp.slice(i-1), triads_cp.slice(i));
        }
        deltaE = 0;
        for (int i=hinge_id-1;i<num_bp-1;i++) {
            deltaE += BPS[i]->eval_delta_energy();
        }

        changed_bps = arma::zeros<arma::ivec>(num_bp-hinge_id);
        int first_id = hinge_id-1;
        for (int i=0;i<num_bp-hinge_id;i++) {
            changed_bps(i) = first_id+i;
        }

        /*
            Contribution due to torque or torsional trap
        */
        deltaE += chain->propose_terminus_rotation(triads->slice(num_bp-1)*Rz, dtheta);

        if (exp(-deltaE) <= choice) {
            chain->set_terminus_rotation(false);
            return false;
        }
        *triads = triads_cp;
        chain->set_terminus_rotation(true);
        return true;
    }
    return true;
}


