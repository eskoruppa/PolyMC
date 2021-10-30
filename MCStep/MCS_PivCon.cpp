#include "MCS_PivCon.h"

/*
    This move never chances the orientation of the first and the last triad.
*/

MCS_PivCon::MCS_PivCon(Chain * ch,const std::vector<long long int> & seedseq, int hingesize_min, int hingesize_max)
: MCStep(ch,seedseq)
{
    move_name="MCS_PivCon";
    /*
        selrange refers to the distance between the the two hinges.
    */
    fac_sigma = 1.2;

    arma::mat covmat = chain->get_avg_cov();
    sigmas     = arma::zeros(3);
    sigmas(0)  = fac_sigma*sqrt(covmat(0,0));
    sigmas(1)  = fac_sigma*sqrt(covmat(1,1));
    sigmas(2)  = fac_sigma*sqrt(covmat(2,2));

    if (hingesize_max >= num_bp/2)   hingesize_max = num_bp/2;
    if (hingesize_min < 2)           hingesize_min = 2;

    int rfrom = 1;
    int rto   = num_bp-1-hingesize_min;

    rge_id1 = 1;
    rge_id2 = num_bp-1;

    decltype(genhinge1.param()) new_range1(rfrom, rto);
    genhinge1.param(new_range1);

    decltype(genhingedist.param()) new_range2(hingesize_min, hingesize_max);
    genhingedist.param(new_range2);

    // intervals of moved segments for exluded volume interactions
    moved_intervals.push_back({0,-1,0});
    moved_intervals.push_back({-1,-1,1000});
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

MCS_PivCon::MCS_PivCon(Chain * ch,const std::vector<long long int> & seedseq, int hingesize_min, int hingesize_max, int range_id1, int range_id2)
: MCStep(ch,seedseq)
/*
    range_id1 and range_id2 specify the 2 pivot points between which the move can change positions and orientations of the chain.
    These monomers themselves can be changed by the move.

    The monomer at range_id1 will never change position, but its attached triad may be rotated.
    The monomer at range_id2 is the last within the move interval whose triad can be rotated.

    All monomers trailing range_id2 will be displaced by the move, but their triads will never be rotated
*/
{
    move_name="MCS_PivCon";
    /*
        selrange refers to the distance between the the two hinges.
    */
    fac_sigma = 1.2;

    arma::mat covmat = chain->get_avg_cov();
    sigmas     = arma::zeros(3);
    sigmas(0)  = fac_sigma*sqrt(covmat(0,0));
    sigmas(1)  = fac_sigma*sqrt(covmat(1,1));
    sigmas(2)  = fac_sigma*sqrt(covmat(2,2));

    if (hingesize_max >= num_bp/2)   hingesize_max = num_bp/2;
    if (hingesize_min < 2)           hingesize_min = 2;

    if (range_id1 < 1) {
        range_id1 = 1;
    }
    range_id2++;
    if (range_id2 >= num_bp) {
        range_id2 = num_bp-1;
    }

    int rfrom = range_id1;
    int rto   = range_id2-hingesize_min;

    rge_id1 = range_id1;
    rge_id2 = range_id2;

    decltype(genhinge1.param()) new_range1(rfrom, rto);
    genhinge1.param(new_range1);

    decltype(genhingedist.param()) new_range2(hingesize_min, hingesize_max);
    genhingedist.param(new_range2);

    // intervals of moved segments for exluded volume interactions
    moved_intervals.push_back({0,-1,0});
    moved_intervals.push_back({-1,-1,1000});
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



MCS_PivCon::~MCS_PivCon() {

}

void MCS_PivCon::update_settings() {
    arma::mat covmat = chain->get_avg_cov();
    sigmas(0)  = fac_sigma*sqrt(covmat(0,0));
    sigmas(1)  = fac_sigma*sqrt(covmat(1,1));
    sigmas(2)  = fac_sigma*sqrt(covmat(2,2));
}


bool MCS_PivCon::MC_move() {

    int    id1,id2,id1m1,id2m1,hinge_dist;
    double deltaE;
    arma::colvec move_Phi;
    arma::mat    rot_mat,Trot_id1,Trot_id2m1;


    id1        = genhinge1(gen);
    hinge_dist = genhingedist(gen);
    id2        = id1 + hinge_dist;

//    if (id2 >= num_bp) {
//        id2 = num_bp-1;
//    }
    if (id2 > rge_id2) {
        id2 = rge_id2;
    }
    id1m1 = id1-1;
    id2m1 = id2-1;

    changed_bps={id1m1,id2m1};

    /*
        Rotate first and last triad of roated segment.
    */

    move_Phi = {normaldist(gen)*sigmas(0),normaldist(gen)*sigmas(1),normaldist(gen)*sigmas(2)};
    rot_mat        = getRotMat(triads->slice(id1m1)*move_Phi);

    Trot_id1   = rot_mat*triads->slice(id1);
    Trot_id2m1 = rot_mat*triads->slice(id2m1);

    BPS[id1m1]->propose_move(triads->slice(id1m1),Trot_id1);
    deltaE  = BPS[id1m1]->eval_delta_energy();

    BPS[id2m1]->propose_move(Trot_id2m1,triads->slice(id2));
    deltaE += BPS[id2m1]->eval_delta_energy();

	if (chain->force_active() == true) {
        arma::colvec trans = pos->col(id1) + rot_mat*(pos->col(id2) - pos->col(id1))  - pos->col(id2);
        deltaE += arma::dot(chain->get_beta_force_vec(),-trans);
	}

	if (exp(-deltaE) <= uniformdist(gen)) {
		return false;
    }
	else {
        triads->slice(id1) = Trot_id1;
        for (int i=id1+1;i<id2m1;i++) {
            triads->slice(i) = rot_mat * triads->slice(i);
            pos->col(i)      = pos->col(i-1) + triads->slice(i-1).col(2)*disc_len;
        }
        triads->slice(id2m1) = Trot_id2m1;
        pos->col(id2m1)      = pos->col(id2m1-1) + triads->slice(id2m1-1).col(2)*disc_len;

        arma::colvec trans = pos->col(id2m1) + triads->slice(id2m1).col(2)*disc_len - pos->col(id2);
        for (int i=id2;i<num_bp;i++) {
            pos->col(i) = pos->col(i) + trans;
        }

        // the first one is permanently set to 0
        moved_intervals[0](EV_TO)   = id1;
        moved_intervals[1](EV_FROM) = id1+1;
        moved_intervals[1](EV_TO)   = id2m1;
        moved_intervals[2](EV_FROM) = id2;
        // the last one is permanently set to num_bp-1
    }
    return true;
}


