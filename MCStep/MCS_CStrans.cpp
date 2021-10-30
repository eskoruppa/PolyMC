#include "MCS_CStrans.h"

/*
    This move never chances the orientation of the first and the last triad.



*/


MCS_CStrans::MCS_CStrans(Chain * ch,const std::vector<long long int> & seedseq, int selrange_min, int selrange_max, arma::colvec dir)
: MCStep(ch,seedseq)
{
    move_name="MCS_CStrs";
    /*
        selrange refers to the distance between the first and the last hinge. (r1 and r3).
        For convenience genhingedist generates the distances between the first and the middle
        hinge, which means that selrange will effectively always be an even number.
    */

    fac = 1.19*2;
    fac = 0.6;
    covmat = chain->get_avg_cov();
    //fac  = 1.76;

    sigma = fac*sqrt(0.5 * (covmat(0,0)+covmat(1,1)) ); //*sqrt(chain->get_T()/300); -> The temperature dependence is already in covmat

    if (arma::norm(dir)>1e-10) {
        fixed_dir=true;
        trans_dir = dir/arma::norm(dir);
    }
    else {
        fixed_dir=false;
    }

    if (selrange_max >= num_bp/2)   selrange_max = num_bp/2;
    if (selrange_min < 2)           selrange_min = 2;

    int hingesize_min = selrange_min/2;
    int hingesize_max = selrange_max/2;

    int rfrom = 1;
    int rto   = num_bp-1-selrange_min;

    decltype(genhinge1.param()) new_range1(rfrom, rto);
    genhinge1.param(new_range1);

    decltype(genhingedist.param()) new_range2(hingesize_min, hingesize_max);
    genhingedist.param(new_range2);

    // intervals of moved segments for exluded volume interactions
    moved_intervals.push_back({0,-1,0});
    moved_intervals.push_back({-1,-1,1000});
    moved_intervals.push_back({-1,-1,1001});
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


MCS_CStrans::~MCS_CStrans() {

}

void MCS_CStrans::update_settings() {
    covmat = chain->get_avg_cov();
    sigma = fac*sqrt(0.5 * (covmat(0,0)+covmat(1,1)) );
}


bool MCS_CStrans::MC_move() {
//    count_step++;

/*
    Declare Variables
*/
    int id1,id2,id3,hingedist;
    int id1n,id2n,id3n;
    arma::colvec r1,r2,r3;
    arma::colvec v1,v2;
    arma::colvec w,d;

    arma::colvec move_Phi;
    arma::mat    rot_mat1,rot_mat2;

    arma::colvec v1p,v2p;
    arma::colvec gamma;

    double dot_gamma_d, discriminant;
    double lambda1,lambda2;

    arma::colvec r2p,r3p;
    arma::colvec r3p1,r3p2,r3_tilde;

    arma::mat    T1,T2m1,T2,T3m1;


/*
    Select Hinges
*/
    hingedist = genhingedist(gen);
    id1 = genhinge1(gen);
    // check if selected segment is within the bounds of the chain
    id3 = id1 + 2*hingedist;
    if  (id3 >= num_bp) {
        hingedist = (num_bp-1-id1)/2;
        id3 = id1 + 2*hingedist;
    }
    id2 = id1 + hingedist;

    id1n = id1-1;
    id2n = id2-1;
    id3n = id3-1;

/*
    Calculate initial position vectors, connecting vectors
    and define the direction of the translation of the tail
    of the chain
*/
    // define initial r's and v's
    r1 = pos->col(id1);
    r2 = pos->col(id2);
    r3 = pos->col(id3);

    v1 = r2-r1;
    v2 = r3-r2;
    w  = r3-r1;

    // define d
    if (fixed_dir) {
        d = trans_dir;
    }
    else {
        d = w/norm(w);
    }


//  This repeat would break the symmetry of the selection probability.
//
//    discriminant = 0;
//    int repeats = 0;
//    while (discriminant <= 0) {
//        if (repeats==50) {
//            cout << "Negative discriminant #" << repeats << " " << discriminant << endl;
//            return false;
//        }
//        repeats++;
//        // Find rotation matrix for the segments between the hinges id1 and id2 (rot_mat1)
//        // The rotation vector is chosen in the frame of the triad id1n. This is not
//        // a problem since the first triad in the chain is never selected as hinge.
//        move_Phi = {normaldist(gen)*sigma,normaldist(gen)*sigma,0};
//        rot_mat1 = getRotMat(triads->slice(id1n)*move_Phi);
//
//        // find new r2 (r2p)
//        v1p = rot_mat1*v1;
//        r2p = r1 + v1p;
//
//    ///////////////////////////////////////////////////////////////
//    // Calculate new r3 (r3p)
//
//    // calculate calculate the 2 possible new positions for r3
//        gamma = r2p-r3;
//        dot_gamma_d  = arma::dot(gamma,d);
//        discriminant = dot_gamma_d*dot_gamma_d+arma::dot(v2,v2)-arma::dot(gamma,gamma);
//    }



    // Find rotation matrix for the segments between the hinges id1 and id2 (rot_mat1)
    // The rotation vector is chosen in the frame of the triad id1n. This is not
    // a problem since the first triad in the chain is never selected as hinge.
    move_Phi = {normaldist(gen)*sigma,normaldist(gen)*sigma,0};
    rot_mat1 = getRotMat(triads->slice(id1n)*move_Phi);

    // find new r2 (r2p)
    v1p = rot_mat1*v1;
    r2p = r1 + v1p;

    ///////////////////////////////////////////////////////////////
    // Calculate new r3 (r3p)

    // calculate calculate the 2 possible new positions for r3
    gamma = r2p-r3;
    dot_gamma_d  = arma::dot(gamma,d);
    discriminant = dot_gamma_d*dot_gamma_d+arma::dot(v2,v2)-arma::dot(gamma,gamma);

    if (discriminant < 0) {
        return false;
    }

    lambda1 = dot_gamma_d + sqrt(discriminant);
    lambda2 = dot_gamma_d - sqrt(discriminant);

    r3p1 = r3 + lambda1*d;
    r3p2 = r3 + lambda2*d;

    // TEST #####################
    assert (std::abs(arma::norm(v2) - arma::norm(r3p1-r2p)) < 1e-12);
    assert (std::abs(arma::norm(v2) - arma::norm(r3p2-r2p)) < 1e-12);
    // ##########################

    // Determine whether random selection is necessary and if not select the r3p that is closest to r3
    r3_tilde = r3-2*arma::dot(v2,d)*d;

    // TEST #####################
    assert (std::abs( arma::norm(r3-r2) - arma::norm(r3_tilde-r2) ) < 1e-12);
    // ##########################

    double d1,d2,d3,d4;
    d1 = arma::norm(r3p1-r3);
    d2 = arma::norm(r3p2-r3);
    d3 = arma::norm(r3p1-r3_tilde);
    d4 = arma::norm(r3p2-r3_tilde);

    bool rand_select=false;
    bool move_to_r3p1 = true;
    if (d1<d2) {
        if (d1>=d3 || d4>=d2 || d4>=d3) {
            rand_select=true;
        }
    }
    else {
        move_to_r3p1 = false;
        // TODO: delete this check. The chance of d1==d2 to occure is negligibly low.
        if (d1==d2) {
            rand_select=true;
            std::cout << "d1 = d2 !!!!!!!!" << std::endl;
        }
        if (d2>=d4 || d3>=d1 || d3>=d4) {
            rand_select=true;
        }
    }

    // if necessary select new position randomly (between the two possibilities)
    if (rand_select) {
        if (uniformdist(gen)<0.5) move_to_r3p1 = true;
        else                      move_to_r3p1 = false;
    }

    if (move_to_r3p1)   r3p = r3p1;
    else                r3p = r3p2;

    ///////////////////////////////////////////////////////////////
    // Find rotation matrix for the segments between the hinges id2 and id3 (rot_mat2)

    v2p = r3p-r2p;

    arma::colvec nv2  = v2/arma::norm(v2);
    arma::colvec nv2p = v2p/arma::norm(v2p);
    arma::colvec axis = arma::cross(nv2,nv2p);
    double theta      = asin(arma::norm(axis));

    /*
    TODO: CHECK THIS THETA!
    */
    if (arma::dot(nv2,nv2p)<0) {
        theta = M_PI-theta; // was 2*M_PI-theta;
    }
//    double theta      = acos(arma::dot(v2/arma::norm(v2),v2p/arma::norm(v2p)));
    axis              = theta*axis/arma::norm(axis);
    rot_mat2          = getRotMat(axis);

//    if (arma::norm(r3p-(r2p+rot_mat2*v2)) > 1e-12) {
//        cout << "difference: "<< arma::norm(r3p-(r2p+rot_mat2*v2)) << endl;
//    }

    // TEST #####################
//    assert( arma::norm(r3p-(r2p+rot_mat2*v2)) < 1e-8 );
    if (arma::norm(r3p-(r2p+rot_mat2*v2)) > 1e-10) {
        std::cout << "Critical Discrepancy!" << std::endl;
    }
    // ##########################


    ///////////////////////////////////////////////////////////////
    // Calculate the energy of the move.

    T1   = rot_mat1*triads->slice(id1);
    T2m1 = rot_mat1*triads->slice(id2n);
    T2   = rot_mat2*triads->slice(id2);
    T3m1 = rot_mat2*triads->slice(id3n);

    BPS[id1n]->propose_move(triads->slice(id1n),T1);
    BPS[id2n]->propose_move(T2m1,T2);
    BPS[id3n]->propose_move(T3m1,triads->slice(id3));

    double deltaE;
    deltaE  = BPS[id1n]->eval_delta_energy();
    deltaE += BPS[id2n]->eval_delta_energy();
    deltaE += BPS[id3n]->eval_delta_energy();

    // Energetic contribution from the force
    if (chain->force_active() == true) {
        deltaE += arma::dot(chain->get_beta_force_vec(),r3-r3p);
	}

	///////////////////////////////////////////////////////////////
	// Assigned changed bps
	changed_bps = {id1n,id2n,id3n};

    ///////////////////////////////////////////////////////////////
    // Metropolis Step

    double choice   = uniformdist(gen);
    if (exp(-deltaE) <= choice) {
//        BPS[id1n]->set_move(false);
//        BPS[id2n]->set_move(false);
//        BPS[id3n]->set_move(false);
		return false;
    }

    ///////////////////////////////////////////////////////////////
    // Move accepted: perform all remaining rotations and translations
//    count_accept++;

    // accept moves
//    BPS[id1n]->set_move(true);
//    BPS[id2n]->set_move(true);
//    BPS[id3n]->set_move(true);

    // Rotate all triads within the move segment
    triads->slice(id1)  = T1;
    triads->slice(id2n) = T2m1;
    triads->slice(id2)  = T2;
    triads->slice(id3n) = T3m1;

    for (int i=id1+1;i<id2n;i++) {
        triads->slice(i) = rot_mat1 * triads->slice(i);
    }
    for (int i=id2+1;i<id3n;i++) {
        triads->slice(i) = rot_mat2 * triads->slice(i);
    }

    // Translate all triads within the move segment
    int im1;
    for (int i=id1+1;i<=id3;i++) {
        im1 = i-1;
        pos->col(i) = pos->col(im1) + triads->slice(im1).col(2)*disc_len;
    }

//    cout << arma::norm(r3p-pos->col(id3)) << endl;
    // TEST #####################
//    assert( arma::norm(r2p-pos->col(id2)) < 1e-10 );
//    assert( arma::norm(r3p-pos->col(id3)) < 1e-10 );
    // ##########################

    // Translate all triads in the tail.
    arma::colvec tail_displacement = pos->col(id3)-r3;
    for (int i=id3+1;i<num_bp;i++) {
        pos->col(i) = pos->col(i)+tail_displacement;
    }

//    if (fixed_dir) {
//        assert( pos->col(num_bp-1)(0) < 0.1 && pos->col(num_bp-1)(1) < 0.1);
//    }

//    cout << "\nMOVED --------------------------------------------" << endl;
//    cout << id1 << " " << id2 << " " << id3 << endl;

    // Intervals of moved segments
    moved_intervals[0](EV_TO)   = id1-1;
    moved_intervals[1](EV_FROM) = id1;
    moved_intervals[1](EV_TO)   = id2;
    moved_intervals[2](EV_FROM) = id2+1;
    moved_intervals[2](EV_TO)   = id3;
    moved_intervals[3](EV_FROM) = id3+1;
    // the last one is permanently set to num_bp-1

    return true;
}


