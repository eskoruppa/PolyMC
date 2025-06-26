#include "MCS_Pivot.h"


MCS_Pivot::MCS_Pivot(Chain * ch,const std::vector<long long int> & seedseq)
: MCStep(ch,seedseq)
{
    move_name="MCS_Pivot";
    // initialize range of random hinge selector
    decltype(genhinge.param()) new_range(0, num_bp-1);
    genhinge.param(new_range);

    fac_theta = 1.2;
    fac_size = 1;

    arma::mat covmat = chain->get_avg_cov();
    sigmas     = arma::zeros(3);
    sigmas(0)  = fac_theta*fac_size*sqrt(covmat(0,0));
    sigmas(1)  = fac_theta*fac_size*sqrt(covmat(1,1));
    sigmas(2)  = fac_theta*fac_size*sqrt(covmat(2,2));

    // intervals of moved segments for exluded volume interactions
    moved_intervals.push_back({0,num_bp-1,0});
    moved_intervals.push_back({-1,num_bp-1,1000});

    suitablility_key[MCS_FORCE_ACTIVE]              = true;  // force_active
    suitablility_key[MCS_TORQUE_ACTIVE]             = false; // torque_active
    suitablility_key[MCS_TOPOLOGY_CLOSED]           = false; // topology_closed
    suitablility_key[MCS_FIXED_LINK]                = false; // link_fixed
    suitablility_key[MCS_FIXED_TERMINI]             = false; // fixed_termini
    suitablility_key[MCS_FIXED_TERMINI_RADIAL]      = false; // fixed_termini_radial
    suitablility_key[MCS_FIXED_FIRST_ORIENTATION]   = false; // fixed_first_orientation
    suitablility_key[MCS_FIXED_LAST_ORIENTATION]    = false; // fixed_last_orientation

    requires_EV_check=true;
}


MCS_Pivot::MCS_Pivot(Chain * ch,const std::vector<long long int> & seedseq, int id1, int id2)
: MCStep(ch,seedseq)
/*
    id1, and id2 specify the (included) boundaries of the interval within which the pivot
    point can be selected. The selected point is the first triad to be rotated. The position
    of the monomer is not changed by the move.
*/
{
    move_name="MCS_Pivot";
    // initialize range of random hinge selector
    if (id2 >= num_bp || id2 < 0) {
        id2 = num_bp - 1;
    }
    decltype(genhinge.param()) new_range(id1,id2);
    genhinge.param(new_range);

    fac_theta = 1.2;
    fac_size = 1;

    arma::mat covmat = chain->get_avg_cov();
    sigmas     = arma::zeros(3);
    sigmas(0)  = fac_theta*fac_size*sqrt(covmat(0,0));
    sigmas(1)  = fac_theta*fac_size*sqrt(covmat(1,1));
    sigmas(2)  = fac_theta*fac_size*sqrt(covmat(2,2));

    // intervals of moved segments for exluded volume interactions
    moved_intervals.push_back({0,num_bp-1,0});
    moved_intervals.push_back({-1,num_bp-1,1000});

    suitablility_key[MCS_FORCE_ACTIVE]              = true;  // force_active
    suitablility_key[MCS_TORQUE_ACTIVE]             = false; // torque_active
    suitablility_key[MCS_TOPOLOGY_CLOSED]           = false; // topology_closed
    suitablility_key[MCS_FIXED_LINK]                = false; // link_fixed
    suitablility_key[MCS_FIXED_TERMINI]             = false; // fixed_termini
    suitablility_key[MCS_FIXED_TERMINI_RADIAL]      = false; // fixed_termini_radial
    suitablility_key[MCS_FIXED_FIRST_ORIENTATION]   = false; // fixed_first_orientation
    suitablility_key[MCS_FIXED_LAST_ORIENTATION]    = false; // fixed_last_orientation

    requires_EV_check=true;
}


MCS_Pivot::~MCS_Pivot() {

}

void MCS_Pivot::update_settings() {
    arma::mat covmat = chain->get_avg_cov();
    sigmas(0)  = fac_theta*fac_size*sqrt(covmat(0,0));
    sigmas(1)  = fac_theta*fac_size*sqrt(covmat(1,1));
    sigmas(2)  = fac_theta*fac_size*sqrt(covmat(2,2));
}

bool MCS_Pivot::MC_move() {
    int    hingeID  = genhinge(gen);
    int    BPS_ID   = hingeID-1;
    double choice   = uniformdist(gen);
    double deltaE   = 0;

    arma::colvec move_Phi = {normaldist(gen)*sigmas(0),normaldist(gen)*sigmas(1),normaldist(gen)*sigmas(2)};
    arma::mat    rot_mat;

    if (hingeID > 0) {
        changed_bps={BPS_ID};

        rot_mat        = getRotMat(triads->slice(hingeID-1)*move_Phi);
        arma::mat Trot = rot_mat*triads->slice(hingeID);

        BPS[BPS_ID]->propose_move(triads->slice(hingeID-1),Trot);
        deltaE = BPS[BPS_ID]->eval_delta_energy();

    }
    else {
        changed_bps.reset();
        rot_mat = getRotMat(move_Phi);
    }

	if (chain->force_active()) {
        // TODO: don_t assume that the first position is at (0,0,0) -> since the first position doesn't change in this move subtracting pos->(0) is technically unnecessary
		arma::colvec e2e_old = pos->col(num_bp-1) - pos->col(0);
		arma::colvec e2e_new = rot_mat*(pos->col(num_bp-1)-pos->col(hingeID)) + pos->col(hingeID) - pos->col(0);
        deltaE += arma::dot(chain->get_beta_force_vec(),e2e_old-e2e_new); //- arma::dot(chain->get_beta_force_vec(),e2e_new);
	}

    if (chain->closure_distance_active()) {
        // std::cout << "closure force active" << std::endl;
        arma::colvec p1 = pos->col(0);
        arma::colvec p2 = pos->col(num_bp-1);
        // old
        deltaE -= chain->eval_closure_distance_betaenergy(p1,p2);
        // new
        p2 = rot_mat*(pos->col(num_bp-1)-pos->col(hingeID)) + pos->col(hingeID);
        deltaE += chain->eval_closure_distance_betaenergy(p1,p2);
    }

    if (chain->closure_angle_active()) {
        // std::cout << "closure angularstiff active" << std::endl;
        arma::colvec tan1 = triads->slice(0).col(2);
        arma::colvec tan2 = triads->slice(num_bp-1).col(2);
        // old
        deltaE -= chain->eval_closure_angle_betaenergy(tan1,tan2);
        // new
        tan2 = rot_mat*tan2;
        deltaE += chain->eval_closure_angle_betaenergy(tan1,tan2);
    }

    if (chain->closure_twist_active()) {
        // std::cout << "closure angularstiff active" << std::endl;
        arma::mat triad1 = triads->slice(num_bp-1);
        arma::mat triad2 = triads->slice(0);
        // old
        deltaE -= chain->eval_closure_twist_betaenergy(triad1,triad2);
        // new
        triad1 = rot_mat*triad1;
        deltaE += chain->eval_closure_twist_betaenergy(triad1,triad2);
    }

	if (exp(-deltaE) <= choice) {
		return false;
    }
	else {
		triads->slice(hingeID) = rot_mat * triads->slice(hingeID);
        for (int j=hingeID+1;j<num_bp;j++) {
            triads->slice(j) = rot_mat * triads->slice(j);
            pos->col(j)      = pos->col(j-1) + triads->slice(j-1).col(2)*disc_len;
        }
        moved_intervals[0](EV_TO) = hingeID;
        moved_intervals[1](EV_FROM) = hingeID+1;
        // the last one is permanently set to num_bp-1
    }
    return true;
}

