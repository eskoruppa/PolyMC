#include "MCS_PivotKink.h"


MCS_PivotKink::MCS_PivotKink(Chain * ch, double A, double D, double B, double x0, double h, double theta, arma::vec T0)
: MCStep(ch), A(A), D(D), B(B), x0(x0), h(h)
{
    move_name="MCS_PivotKink";
    // initialize range of random hinge selector
    decltype(genhinge.param()) new_range(0, num_bp-1);
    genhinge.param(new_range);

//    sigfac = 4.130;
    fac_theta = 1.2;
    fac_size = 1;

    xcon     = x0*0.5*(1+2*h/(A*x0*x0));
    costheta = std::cos(theta);
    sintheta = std::sin(theta);

    A = A*0.5/disc_len;
    D = D*0.5/disc_len;
    B = B*0.5/disc_len;

    R0   = getRotMat(T0);
    R0_T = R0.t();

    cout << R0 << endl;

//    arma::mat covmat = chain->get_avg_cov();

    double fac = 0.25;

    sigmas     = arma::zeros(3);
    sigmas(0)  = 1*fac;
    sigmas(1)  = 1*fac;
    sigmas(2)  = 1*fac;

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

MCS_PivotKink::~MCS_PivotKink() {

}

void MCS_PivotKink::update_settings() {
    arma::mat covmat = chain->get_avg_cov();
    sigmas(0)  = fac_theta*fac_size*sqrt(covmat(0,0));
    sigmas(1)  = fac_theta*fac_size*sqrt(covmat(1,1));
    sigmas(2)  = fac_theta*fac_size*sqrt(covmat(2,2));
}

bool MCS_PivotKink::MC_move() {

    int    hingeID  = genhinge(gen);
    int    BPS_ID   = hingeID-1;
    double choice   = uniformdist(gen);
    double deltaE=0;

    arma::colvec move_Phi = {normaldist(gen)*sigmas(0),normaldist(gen)*sigmas(1),normaldist(gen)*sigmas(2)};
    arma::mat    rot_mat;

    if (hingeID > 0) {
        changed_bps={BPS_ID};

        rot_mat        = getRotMat(triads->slice(hingeID-1)*move_Phi);
//        rot_mat        = getRotMat(ExtractTheta( getRotMat(move_Phi)*triads->slice(hingeID-1).t()));
//        rot_mat        = getRotMat(ExtractTheta(triads->slice(hingeID-1)*getRotMat(move_Phi)));
        arma::mat Trot = rot_mat*triads->slice(hingeID);

//        BPS[BPS_ID]->propose_move(triads->slice(hingeID-1),Trot);
//        deltaE = BPS[BPS_ID]->eval_delta_energy();

        arma::colvec theta_old = ExtractTheta( R0_T * triads->slice(hingeID-1).t() * triads->slice(hingeID) );
//        cout << theta_old.t();
//        cal_energy(theta_old);

        arma::colvec theta_new = ExtractTheta( R0_T * triads->slice(hingeID-1).t() * Trot );
        deltaE = cal_energy(theta_new) - cal_energy(theta_old);
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
//        if (hingeID > 0){
//            BPS[BPS_ID]->set_move(false);
//        }
		return false;
    }
	else {
//        count_accept++;
		triads->slice(hingeID) = rot_mat * triads->slice(hingeID);
        for (int j=hingeID+1;j<num_bp;j++) {
            triads->slice(j) = rot_mat * triads->slice(j);
            pos->col(j)      = pos->col(j-1) + triads->slice(j-1).col(2)*disc_len;
        }

//		if (hingeID > 0) {
//            BPS[BPS_ID]->set_move(true);
//		}

        moved_intervals[0](EV_TO) = hingeID;
        moved_intervals[1](EV_FROM) = hingeID+1;
        // the last one is permanently set to num_bp-1
    }
    return true;
}

double MCS_PivotKink::cal_energy(arma::colvec Theta) {
    double x,y,z;
    x = Theta(1)*costheta-Theta(2)*sintheta;
    y = Theta(1)*sintheta+Theta(2)*costheta;
    z = Theta(0);

//    cout << "x = " << x << " (x_con = " << xcon << ")" << endl;
//    std::exit(1);

    double E = D*y*y + B*z*z;
    if (x>xcon) {
        E += A*x*x;
    }
    else {
//        cout << "############" << endl;
//        cout << "x = " << x << endl;
//        cout << Theta.t();
//        cout << "costheta = " << costheta << endl;
//        cout << "sintheta = " << sintheta << endl;
        double dx = x-x0;
        E += A*dx*dx+h;
    }
    return E;
}


