#include "MCS_Twist.h"


MCS_Twist::MCS_Twist(Chain * ch,const std::vector<long long int> & seedseq)
: MCStep(ch,seedseq)
{
    move_name = "MCS_Twist";

    sigfac = 2;
    arma::mat covmat = chain->get_avg_cov();
    sigma = sigfac*sqrt(covmat(2,2));

    // initialize range of random hinge selector
    int sf=0;
    int sl=num_bp-1;
//    if (!chain->topology_closed()) {
//        sf=1;
//        sl=num_bp-2;
//    }

    if (chain->fixed_first_orientation()) sf=1;
    if (chain->fixed_last_orientation())  sl=num_bp-2;
    decltype(genhinge.param()) new_range(sf, sl);
    genhinge.param(new_range);

    suitablility_key[MCS_FORCE_ACTIVE]              = true;  // force_active
    suitablility_key[MCS_TORQUE_ACTIVE]             = true;  // torque_active
    suitablility_key[MCS_TOPOLOGY_CLOSED]           = true;  // topology_closed
    suitablility_key[MCS_FIXED_LINK]                = true;  // link_fixed
    suitablility_key[MCS_FIXED_TERMINI]             = true;  // fixed_termini
    suitablility_key[MCS_FIXED_TERMINI_RADIAL]      = true;  // fixed_termini_radial
    suitablility_key[MCS_FIXED_FIRST_ORIENTATION]   = true;  // fixed_first_orientation
    suitablility_key[MCS_FIXED_LAST_ORIENTATION]    = true;  // fixed_last_orientation

    requires_EV_check=false;
}


MCS_Twist::MCS_Twist(Chain * ch,const std::vector<long long int> & seedseq, int range_id1, int range_id2)
: MCStep(ch,seedseq)
/*
    range_id1 and range_id2 define the (included) boundary of the range in which triads
    can be rotated by the move.
    For closed chains the range_id2 can be set to a smaller value than range_id1, which
    defines an interval accross the periodic boundary.
*/
{
    move_name = "MCS_Twist";

    sigfac = 2;
    arma::mat covmat = chain->get_avg_cov();
    sigma = sigfac*sqrt(covmat(2,2));

    // initialize range of random hinge selector

    range_id1 = pmod(range_id1,num_bp);
    range_id2 = pmod(range_id2,num_bp);

    if (range_id2 < range_id1) {
        if (chain->topology_closed()) {
            range_id2 = range_id2 + num_bp;
        }
        else {
            std::cout << "Error invalid Range selection in MCS_Twist. For open chains range_id2 needs to be larger or equal to range_id1!" << std::endl;
            std::exit(0);
        }
    }

    if (chain->fixed_first_orientation() && range_id1==0       ) range_id1=1;
    if (chain->fixed_last_orientation()  && range_id2==num_bp-1) range_id2=num_bp-2;
    decltype(genhinge.param()) new_range(range_id1, range_id2);
    genhinge.param(new_range);

    suitablility_key[MCS_FORCE_ACTIVE]              = true;  // force_active
    suitablility_key[MCS_TORQUE_ACTIVE]             = true;  // torque_active
    suitablility_key[MCS_TOPOLOGY_CLOSED]           = true;  // topology_closed
    suitablility_key[MCS_FIXED_LINK]                = true;  // link_fixed
    suitablility_key[MCS_FIXED_TERMINI]             = true;  // fixed_termini
    suitablility_key[MCS_FIXED_TERMINI_RADIAL]      = true;  // fixed_termini_radial
    suitablility_key[MCS_FIXED_FIRST_ORIENTATION]   = true;  // fixed_first_orientation
    suitablility_key[MCS_FIXED_LAST_ORIENTATION]    = true;  // fixed_last_orientation

    requires_EV_check=false;
}


MCS_Twist::~MCS_Twist() {

}

void MCS_Twist::update_settings() {
    arma::mat covmat = chain->get_avg_cov();
    sigma = sigfac*sqrt(covmat(2,2));
}

bool MCS_Twist::MC_move() {

    int    id       = pmod(genhinge(gen),num_bp);
    double choice   = uniformdist(gen);
    double deltaE;

    int idm1 = pmod(id-1,num_bp);
    int idp1 = pmod(id+1,num_bp);

    arma::mat Trot   = triads->slice(id)*Rotz(normaldist(gen)*sigma);

    BPS[idm1]->propose_move(triads->slice(idm1), Trot);
    BPS[id  ]->propose_move(Trot, triads->slice(idp1));

    deltaE  = BPS[idm1]->eval_delta_energy();
    deltaE += BPS[id  ]->eval_delta_energy();


    if (exp(-deltaE) <= choice) {
		return false;
    }
    triads->slice(id) = Trot;
    changed_bps = {idm1,id};
    return true;
}


