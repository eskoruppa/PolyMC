#include "MCS_Pivot2d.h"


MCS_Pivot2d::MCS_Pivot2d(Chain * ch,const std::vector<long long int> & seedseq, arma::colvec normal_to_plane)
: MCStep(ch,seedseq),normal(normal_to_plane)
{
    move_name="MCS_Pivot2d";
    // initialize range of random hinge selector
    decltype(genhinge.param()) new_range(0, num_bp-1);
//    decltype(genhinge.param()) new_range(num_bp-5, num_bp-1);
    genhinge.param(new_range);

//    sigfac = 4.130;
    fac    = 5;
    arma::mat covmat = chain->get_avg_cov();
    sigma  = sqrt(sqrt(covmat(0,0)*covmat(0,0)+covmat(1,1)*covmat(1,1))) *  fac*sqrt(chain->get_T()/300)*std::sqrt(disc_len/0.34);

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

MCS_Pivot2d::MCS_Pivot2d(Chain * ch,const std::vector<long long int> & seedseq, arma::colvec normal_to_plane, int id1, int id2)
: MCStep(ch,seedseq),normal(normal_to_plane)
/*
    id1, and id2 specify the (included) boundaries of the interval within which the pivot
    point can be selected. The selected point is the first triad to be rotated. The position
    of the monomer is not changed by the move.
*/
{
    move_name="MCS_Pivot2d";
    // initialize range of random hinge selector
    if (id2 >= num_bp || id2 < 0) {
        id2 = num_bp - 1;
    }
    decltype(genhinge.param()) new_range(id1,id2);
    genhinge.param(new_range);

//    sigfac = 4.130;
    fac    = 5;
    arma::mat covmat = chain->get_avg_cov();
    sigma  = sqrt(sqrt(covmat(0,0)*covmat(0,0)+covmat(1,1)*covmat(1,1))) *  fac*sqrt(chain->get_T()/300)*std::sqrt(disc_len/0.34);

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

MCS_Pivot2d::~MCS_Pivot2d() {

}

void MCS_Pivot2d::update_settings() {
    arma::mat covmat = chain->get_avg_cov();
    sigma  = sqrt(sqrt(covmat(0,0)*covmat(0,0)+covmat(1,1)*covmat(1,1))) *  fac*sqrt(chain->get_T()/300)*std::sqrt(disc_len/0.34);
}

bool MCS_Pivot2d::MC_move() {

    int    hingeID  = genhinge(gen);
    int    BPS_ID   = hingeID-1;
    double choice   = uniformdist(gen);
    double deltaE=0;

    arma::colvec move_Phi = normal*sigma*normaldist(gen);
    arma::mat    rot_mat;

    if (hingeID > 0) {
        changed_bps={BPS_ID};

        rot_mat        = getRotMat(move_Phi);
//        rot_mat        = getRotMat(ExtractTheta( getRotMat(move_Phi)*triads->slice(hingeID-1).t()));
//        rot_mat        = getRotMat(ExtractTheta(triads->slice(hingeID-1)*getRotMat(move_Phi)));
        arma::mat Trot = rot_mat*triads->slice(hingeID);

        BPS[BPS_ID]->propose_move(triads->slice(hingeID-1),Trot);
        deltaE = BPS[BPS_ID]->eval_delta_energy();

    }
    else {
        changed_bps.reset();
        rot_mat = getRotMat(move_Phi);
    }
	if (chain->force_active() == true) {
        // TODO: don_t assume that the first position is at (0,0,0)
		arma::colvec ETE_old = pos->col(num_bp-1);
		arma::colvec ETE_new = rot_mat*(pos->col(num_bp-1)-pos->col(hingeID)) + pos->col(hingeID);
        deltaE += arma::dot(chain->get_beta_force_vec(),ETE_old-ETE_new); //- arma::dot(chain->get_beta_force_vec(),ETE_new);
	}

	if (exp(-deltaE) <= choice) {
		return false;
    }
	else {
//        count_accept++;
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


