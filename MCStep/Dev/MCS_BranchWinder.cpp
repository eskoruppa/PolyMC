#include "MCS_BranchWinder.h"


MCS_BranchWinder::MCS_BranchWinder(Chain * ch, double seg_size, double max_dist, double min_contour_dist)
: MCStep(ch),
max_dist(max_dist)
{
    move_name = "MCS_BrnWd";
    /*
    Initializes the controling variables
    */
    double fac  = 0.18;
    fac = 0.15;
    fac = 0.25;
    sigma = fac*std::sqrt(disc_len/0.34); //sqrt(chain->get_T()/300)*;

    // amount of beads in the rotated segment
    n_beads_in_seg = seg_size/disc_len;
    if (n_beads_in_seg <= 0) n_beads_in_seg = 1;

    cout << "\n";
    cout << "BranchWinder segment size set to " << n_beads_in_seg << endl;

    // minimum bead distance between the two hinges
    min_contour_bead_dist = min_contour_dist/disc_len+2*n_beads_in_seg;
    if (min_contour_bead_dist<MCSBW_MIN_HINGE_DIST) min_contour_bead_dist=MCSBW_MIN_HINGE_DIST;



    /*
        Set the range of the random generator selecting hinge A.
    */
    decltype(genhingeA.param()) new_range(0, num_bp-1);
    genhingeA.param(new_range);

    closed = chain->topology_closed();

    /*
        intervals of moved segments for exluded volume interactions
        Typically this will contain 5 intervals, unless one of the
        intervals starts at the first bead or terminates at the last bead.
    */
    moved_intervals.push_back({-1,-1,-1});
    moved_intervals.push_back({-1,-1,-1});
    moved_intervals.push_back({-1,-1,-1});
    moved_intervals.push_back({-1,-1,-1});

    changed_bps = arma::zeros<arma::ivec>(3);

    suitablility_key[MCS_FORCE_ACTIVE]              = true;  // force_active
    suitablility_key[MCS_TORQUE_ACTIVE]             = true;  // torque_active
    suitablility_key[MCS_TOPOLOGY_CLOSED]           = true;  // topology_closed
    suitablility_key[MCS_FIXED_LINK]                = true;  // link_fixed
    suitablility_key[MCS_FIXED_TERMINI]             = true;  // fixed_termini
    suitablility_key[MCS_FIXED_TERMINI_RADIAL]      = true;  // fixed_termini_radial
    suitablility_key[MCS_FIXED_FIRST_ORIENTATION]   = true;  // fixed_first_orientation
    suitablility_key[MCS_FIXED_LAST_ORIENTATION]    = true;  // fixed_last_orientation

//    requires_EV_check=false;
}


MCS_BranchWinder::~MCS_BranchWinder() {

}

void MCS_BranchWinder::update_settings() {

}


bool MCS_BranchWinder::MC_move() {

    int idA,idB,idC,idD;
    bool moveAC;

    int hingefind_retries;

    /*
        Find the hinges A and B
    */
    idB = -1;
    hingefind_retries = 0;
    while (idB==-1 && hingefind_retries < MCSBW_HINGEFIND_RETRIES) {
        idA = genhingeA(gen);

        #if MCSBW_CHOOSE_CLOSEST == 1
        idB = get_partner_closest(idA);
        #else
        idB = get_partner_random(idA);
        #endif

        hingefind_retries++;
    }
    if (idB==-1) {
        return false;
    }

//    if (get_partner_closest(idA)==-1) {
//        cout << "Closest partner algorithm finds nothing!" << endl;
//        std::exit(1);
//    }
//    if (get_partner_closest(idA) != get_partner_closest_test(idA)) {
//        cout << "Problem with finding closest bead!" << endl;
//        std::exit(1);
//    }


    /*
        For ease of calculation assure that the order is always
        first idA and then idB
    */
    if (!direction_is_forward(idA,idB)) {
        int temp = idA;
        idA = idB;
        idB = temp;
    }

//    cout << "Forwards dist = " << pmod(idB-idA,num_bp) << endl;
//    cout << "Backward dist = " << pmod(idA-idB,num_bp) << endl;

    /*
        Calculate Hinges C and D
    */
    idC = pmod(idA+n_beads_in_seg,num_bp);
    idD = pmod(idB-n_beads_in_seg,num_bp);

    /*
        Determine whether the vector AC or the BD is moved.
    */
    if (uniformdist(gen)<0.5)   moveAC = true;
    else                        moveAC = false;

//    moveAC = true;

    arma::colvec rA,rB,rC,rD;
    if (moveAC) {
        rA = pos->col(idB);
        rB = pos->col(idA);
        rC = pos->col(idD);
        rD = pos->col(idC);
    }
    else {
        rA = pos->col(idA);
        rB = pos->col(idB);
        rC = pos->col(idC);
        rD = pos->col(idD);
    }

    /*
        Calculate the rotation
    */

    arma::colvec vec_v1,vec_u, uvec_n;
//    arma::colvec vec_v2;
    arma::colvec rM,rDp;
    double theta;

    vec_v1 = rD-rB;
//    vec_v2 = rC-rA;
    uvec_n = rC-rB;
    uvec_n = uvec_n/arma::norm(uvec_n);

    rM    = rB + arma::dot(uvec_n,vec_v1)*uvec_n;
    vec_u = rD-rM;

    theta  = normaldist(gen)*sigma;
    rDp = rM + getRotMat(theta*uvec_n)*vec_u;


///////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////
/////// TEST //////////////////////////////////////////////
//    double distCD  = arma::norm(rD-rC);
//    double distCDp = arma::norm(rDp-rC);
//    double distBD  = arma::norm(rD-rB);
//    double distBDp = arma::norm(rDp-rB);
//
//    if (distCD - distCDp >1e-12 || distBD - distBDp >1e-12) {
//        cout << "Error in BranchWinder rotation" << endl;
//        cout << arma::dot(vec_u,uvec_n) << endl;
//        cout << "delta CD: " << distCD - distCDp << endl;
//        cout << "delta BD: " << distBD - distBDp << endl;
//    }

///////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////

    /*
        Calculate reference triads
    */

    arma::colvec x,y,z;
    arma::colvec ref_pos1, ref_pos2;
    arma::mat Tref_1(3,3);
    arma::mat Tref_2(3,3);
    arma::mat Rb;

    ref_pos1 = 0.5*(rC+rD);
    ref_pos2 = 0.5*(rC+rDp);

    y  = rD-rC;
    y  = y / arma::norm(y);
    z  = ref_pos1 - 0.5*(rA+rB);
    z  = z - arma::dot(z,y)*y;
    z  = z/arma::norm(z);
    x  = arma::cross(y,z);

    Tref_1.col(0) = x;
    Tref_1.col(1) = y;
    Tref_1.col(2) = z;

    y  = rDp-rC;
    y  = y / arma::norm(y);
    z  = ref_pos2 - 0.5*(rA+rB);
    z  = z - arma::dot(z,y)*y;
    z  = z/arma::norm(z);
    x = arma::cross(y,z);

    Tref_2.col(0) = x;
    Tref_2.col(1) = y;
    Tref_2.col(2) = z;

    Rb = Tref_2*Tref_1.t();



//    cout << "theta = " << theta*180/M_PI << endl;
//    arma::colvec Omega = ExtractTheta(R);
//    cout << Omega.t()*180/M_PI;
//    if (Omega.has_nan()) {
//        cout << "Nan detected" << endl;
//        cout << Omega.t()*180/M_PI;
//        cout << R << endl;
//        cout << Tref_1 << endl;
//        cout << Tref_2 << endl;
//        cout << (rD-rC).t();
//        cout << idA << " - " << idB << endl;
//        cout << idC << " - " << idD << endl;
//    }
//    arma::colvec Cconnect, Dconnect;
//    Cconnect = rC - ref_pos1;
//    Dconnect = rD - ref_pos1;
//    arma::colvec Cp_test, Dp_test;
//    Cp_test = ref_pos2+R*Cconnect;
//    Dp_test = ref_pos2+R*Dconnect;
//
//    cout << "Rotation test " << endl;
//    cout << rC-Cp_test << endl;
//    cout << rDp-Dp_test << endl;
//
//    cout << "Rotation test " << endl;
//    cout << rDp-rB+R*vec_v1 << endl;

    /*
        Find the rotation matrix for the segment v1
    */

    arma::mat Rv1;
    arma::colvec Rv1_rotax;
    double len_Rv1_rotax;
    double Rv1_phi;
    arma::colvec vec_v1p, uvec_v1, uvec_v1p;
    double len_v1;

    len_v1  = arma::norm(vec_v1);
    vec_v1p = rDp-rB;
    uvec_v1  = vec_v1/len_v1;
    uvec_v1p = vec_v1p/len_v1;

    Rv1_rotax = arma::cross(uvec_v1,uvec_v1p);
    len_Rv1_rotax = arma::norm(Rv1_rotax);
    Rv1_phi = std::asin(len_Rv1_rotax);
    if (arma::dot(uvec_v1,uvec_v1p)<0) {
        Rv1_phi = M_PI-Rv1_phi;
    }
    Rv1 = getRotMat(Rv1_phi*Rv1_rotax/len_Rv1_rotax);


//    cout << "Rotation test " << endl;
//    cout << (rDp-rB-Rv1*vec_v1).t() << endl;


    /*
        Evaluate the elastic Energy of the suggested configuration
    */

    double deltaE = 0;

    if (moveAC) {
        int idAm1 = pmod(idA-1,num_bps);
        int idCm1 = pmod(idC-1,num_bps);
        int idDm1 = pmod(idD-1,num_bps);

        BPS[idAm1]-> propose_move( triads->slice(idAm1)     , Rv1*triads->slice(idA)   );
        BPS[idCm1]-> propose_move( Rv1*triads->slice(idCm1) , Rb*triads->slice(idC));
        BPS[idDm1]-> propose_move(Rb*triads->slice(idDm1)   , triads->slice(idD));

        deltaE += BPS[idAm1]->eval_delta_energy();
        deltaE += BPS[idCm1]->eval_delta_energy();
        deltaE += BPS[idDm1]->eval_delta_energy();

        changed_bps(0) = idAm1;
        changed_bps(1) = idCm1;
        changed_bps(2) = idDm1;
    }
    else {
        int idCm1 = pmod(idC-1,num_bps);
        int idDm1 = pmod(idD-1,num_bps);
        int idBm1 = pmod(idB-1,num_bps);

        BPS[idCm1]-> propose_move(    triads->slice(idCm1) , Rb *triads->slice(idC)   );
        BPS[idDm1]-> propose_move(Rb* triads->slice(idDm1) , Rv1*triads->slice(idD)   );
        BPS[idBm1]-> propose_move(Rv1*triads->slice(idBm1) , triads->slice(idB) );

        deltaE += BPS[idCm1]->eval_delta_energy();
        deltaE += BPS[idDm1]->eval_delta_energy();
        deltaE += BPS[idBm1]->eval_delta_energy();

        changed_bps(0) = idCm1;
        changed_bps(1) = idDm1;
        changed_bps(2) = idBm1;
    }

    if (exp(-deltaE) <= uniformdist(gen)) {
		return false;
    }
    // If Move is accepted rotate full segments
	else {

//        double E_old = chain->recal_energy();


        if (moveAC) {

//            cout << "######################################################" << endl;
//            cout << "MoveAC" << endl;

            int id,idp1;
            int num;
            idp1 = idA;

            num = pmod(idC-idA,num_bp);
            for (int i=0;i<num;i++) {
                id = idp1;
                idp1 = pmod(id+1,num_bp);
                triads->slice(id) = Rv1*triads->slice(id);
                pos->col(idp1)    = pos->col(id) + triads->slice(id).col(2)*disc_len;
            }

//            cout << "Test pos C " << endl;
//            cout << (pos->col(idC) - rDp).t() << endl;

            num = pmod(idD-idC,num_bp);
            for (int i=0;i<num;i++) {
                id = idp1;
                idp1 = pmod(id+1,num_bp);
                triads->slice(id) = Rb*triads->slice(id);
                pos->col(idp1)    = pos->col(id) + triads->slice(id).col(2)*disc_len;
            }

//            cout << "Test pos D " << endl;
//            cout << pos->col(idD) - rC << endl;


        }
        else {

//            cout << "######################################################" << endl;
//            cout << "MoveBD" << endl;

            int id,idp1;
            int num;
            idp1 = idC;

            num = pmod(idD-idC,num_bp);
            for (int i=0;i<num;i++) {
                id = idp1;
                idp1 = pmod(id+1,num_bp);
                triads->slice(id) = Rb*triads->slice(id);
                pos->col(idp1)    = pos->col(id) + triads->slice(id).col(2)*disc_len;
            }

//            cout << "Test pos D " << endl;
//            cout << (pos->col(idD) - rDp).t() << endl;

            num = pmod(idB-idD,num_bp);
            for (int i=0;i<num;i++) {
                id = idp1;
                idp1 = pmod(id+1,num_bp);
                triads->slice(id) = Rv1*triads->slice(id);
                pos->col(idp1)    = pos->col(id) + triads->slice(id).col(2)*disc_len;
            }

//            cout << "Test pos B " << endl;
//            cout << pos->col(idB) - rB << endl;


        }

        /*
            CHECK IF I AM CALUCLATING THE RIGHT DELTA E BY CALCULATING THE DIFFERENCE IN ENERGY OF THE ENTIRE CHAIN
        */
//        double E_new = chain->recal_energy();
//        if (std::abs(deltaE-E_new+E_old)>1e-12) {
//            cout << "delta E = " << deltaE << " - " << E_new-E_old << endl;
//            cout << "Energy discrepancy: " << std::abs(deltaE-E_new+E_old) << endl;
//            assert(std::abs(deltaE-E_new+E_old)<1e-12);
//        }
//        cout << "Energy discrepancy: " << std::abs(deltaE-E_new+E_old) << endl;
    }

    /*
        Define the Moved Intervals
    */

    if (moveAC) {
        /*
            Moved segment AC
        */

        if (idA < idC && idA < idD) {
            moved_intervals[0](0) = 0;
            moved_intervals[0](1) = idA;
            moved_intervals[0](2) = 0;

            moved_intervals[1](0) = idA+1;
            moved_intervals[1](1) = idC;
            moved_intervals[1](2) = 1000;

            moved_intervals[2](0) = idC+1;
            moved_intervals[2](1) = idD-1;
            moved_intervals[2](2) = 1001;

            moved_intervals[3](0) = idD;
            moved_intervals[3](1) = num_bp-1;
            moved_intervals[3](2) = 0;
        }
        else {
            if (idC < idA) {
                moved_intervals[0](0) = 0;
                moved_intervals[0](1) = idC;
                moved_intervals[0](2) = 1000;

                moved_intervals[1](0) = idC+1;
                moved_intervals[1](1) = idD-1;
                moved_intervals[1](2) = 1001;

                moved_intervals[2](0) = idD;
                moved_intervals[2](1) = idA;
                moved_intervals[2](2) = 0;

                moved_intervals[3](0) = idA+1;
                moved_intervals[3](1) = num_bp-1;
                moved_intervals[3](2) = 1000;
            }
            else {
                moved_intervals[0](0) = 0;
                moved_intervals[0](1) = idD-1;
                moved_intervals[0](2) = 1001;

                moved_intervals[1](0) = idD;
                moved_intervals[1](1) = idA;
                moved_intervals[1](2) = 0;

                moved_intervals[2](0) = idA+1;
                moved_intervals[2](1) = idC;
                moved_intervals[2](2) = 1000;

                moved_intervals[3](0) = idC+1;
                moved_intervals[3](1) = num_bp-1;
                moved_intervals[3](2) = 1001;
            }
        }
    }
    else {
        /*
            Moved segment BD
        */
        if (idC < idD && idC < idB) {
            moved_intervals[0](0) = 0;
            moved_intervals[0](1) = idC;
            moved_intervals[0](2) = 0;

            moved_intervals[1](0) = idC+1;
            moved_intervals[1](1) = idD-1;
            moved_intervals[1](2) = 1001;

            moved_intervals[2](0) = idD;
            moved_intervals[2](1) = idB-1;
            moved_intervals[2](2) = 1000;

            moved_intervals[3](0) = idB;
            moved_intervals[3](1) = num_bp-1;
            moved_intervals[3](2) = 0;
        }
        else {
            if (idD < idC) {
                moved_intervals[0](0) = 0;
                moved_intervals[0](1) = idD-1;
                moved_intervals[0](2) = 1001;

                moved_intervals[1](0) = idD;
                moved_intervals[1](1) = idB-1;
                moved_intervals[1](2) = 1000;

                moved_intervals[2](0) = idB;
                moved_intervals[2](1) = idC;
                moved_intervals[2](2) = 0;

                moved_intervals[3](0) = idC+1;
                moved_intervals[3](1) = num_bp-1;
                moved_intervals[3](2) = 1001;
            }
            else {
                moved_intervals[0](0) = 0;
                moved_intervals[0](1) = idB-1;
                moved_intervals[0](2) = 1000;

                moved_intervals[1](0) = idB;
                moved_intervals[1](1) = idC;
                moved_intervals[1](2) = 0;

                moved_intervals[2](0) = idC+1;
                moved_intervals[2](1) = idD-1;
                moved_intervals[2](2) = 1001;

                moved_intervals[3](0) = idD;
                moved_intervals[3](1) = num_bp-1;
                moved_intervals[3](2) = 1000;
            }
        }
    }
    return true;
}


int MCS_BranchWinder::get_partner_random(int idA) {

/*
    This can be made much more efficient by using Proximity Clusters
*/
    arma::colvec rA,rB;
    double dist;

    rA = pos->col(idA);
    vector<int> candidates;

    int id_lower = 0;
    int id_upper = num_bp-1;

    if (closed) {
        if (num_bp-1-idA<min_contour_bead_dist) {
            id_lower = pmod(idA+min_contour_bead_dist,num_bp);
        }
        if (idA<min_contour_bead_dist) {
            id_upper = pmod(idA-min_contour_bead_dist,num_bp);
        }
    }

    int idB=id_lower;
    int delta_id;

    while (idB <= idA-min_contour_bead_dist) {

        rB = pos->col(idB);
        dist = arma::norm(rA-rB);
        if (dist<max_dist) {
            candidates.push_back(idB);
        }
//        idB += dist/disc_len;
//        idB++;
        delta_id = (dist-max_dist)/dist;
        if (delta_id<1) delta_id = 1;
        idB += delta_id;

//        if (dist < disc_len) {
//            cout << "Beads  too close (first loop)" << endl;
//            cout << "idA  = " << idA << endl;
//            cout << "idB  = " << idB << endl;
//            cout << "dist = " << dist << endl;
//            cout << "dlen = " << disc_len << endl;
//
//            if (EV->check_overlap()) {
//                cout << "EV overlap found!" << endl;
//            }
//            else {
//                cout << "No EV overlap detected!" << endl;
//            }
////            return -1;
//            std::exit(1);
//        }
    }

    idB=idA+min_contour_bead_dist;
    while (idB < id_upper) {
        rB = pos->col(idB);
        dist = arma::norm(rA-rB);
        if (dist<max_dist) {
            candidates.push_back(idB);
        }
//        idB += dist/disc_len;
//        idB++;
        delta_id = (dist-max_dist)/dist;
        if (delta_id<1) delta_id = 1;
        idB += delta_id;

//        if (dist < disc_len) {
//            cout << "Beads  too close (second loop)" << endl;
//            cout << "idA  = " << idA << endl;
//            cout << "idB  = " << idB << endl;
//            cout << "dist = " << dist << endl;
//            cout << "dlen = " << disc_len << endl;
////            return -1;
//            std::exit(1);
//        }
    }

    decltype(selecthingeB.param()) new_range(0, candidates.size()-1);
    selecthingeB.param(new_range);

    if (candidates.size()==0) {
        return -1;
    }
    return candidates[selecthingeB(gen)];
}

int MCS_BranchWinder::get_partner_closest(int idA) {

/*
    This can be made much more efficient by using Proximity Clusters
*/
    arma::colvec rA,rB;
    double dist;

    rA = pos->col(idA);
    int id_lower = 0;
    int id_upper = num_bp-1;

    if (closed) {
        if (num_bp-1-idA<min_contour_bead_dist) {
            id_lower = pmod(idA+min_contour_bead_dist,num_bp);
        }
        if (idA<min_contour_bead_dist) {
            id_upper = pmod(idA-min_contour_bead_dist,num_bp);
        }
    }

    double closest_dist = max_dist;
    int    closest_id   = -1;

    int idB=id_lower;
    int delta_id;

    while (idB <= idA-min_contour_bead_dist) {

        rB = pos->col(idB);
        dist = arma::norm(rA-rB);
        if (dist<closest_dist) {
            closest_id = idB;
            closest_dist = dist;
        }
//        idB += dist/disc_len;
//        idB++;
        delta_id = (dist-max_dist)/dist;
        if (delta_id<1) delta_id = 1;
        idB += delta_id;
    }

    idB=idA+min_contour_bead_dist;
    while (idB < id_upper) {
        rB = pos->col(idB);
        dist = arma::norm(rA-rB);
        if (dist<closest_dist) {
            closest_id = idB;
            closest_dist = dist;
        }
//        idB += dist/disc_len;
//        idB++;
        delta_id = (dist-max_dist)/dist;
        if (delta_id<1) delta_id = 1;
        idB += delta_id;

    }
    return closest_id;
}

int MCS_BranchWinder::get_partner_closest_test(int idA) {

/*
    This can be made much more efficient by using Proximity Clusters
*/
    arma::colvec rA,rB;
    double dist;

    rA = pos->col(idA);
    int id_lower = 0;
    int id_upper = num_bp-1;

    if (closed) {
        if (num_bp-1-idA<min_contour_bead_dist) {
            id_lower = pmod(idA+min_contour_bead_dist,num_bp);
        }
        if (idA<min_contour_bead_dist) {
            id_upper = pmod(idA-min_contour_bead_dist,num_bp);
        }
    }

    double closest_dist = max_dist;
    int    closest_id   = -1;

    int idB=id_lower;
    int delta_id;

    while (idB <= idA-min_contour_bead_dist) {

        rB = pos->col(idB);
        dist = arma::norm(rA-rB);
        if (dist<closest_dist) {
            closest_id = idB;
            closest_dist = dist;
        }
//        idB += dist/disc_len;
        idB++;
//        delta_id = (dist-max_dist)/dist;
//        if (delta_id<1) delta_id = 1;
//        idB += delta_id;
    }

    idB=idA+min_contour_bead_dist;
    while (idB < id_upper) {
        rB = pos->col(idB);
        dist = arma::norm(rA-rB);
        if (dist<closest_dist) {
            closest_id = idB;
            closest_dist = dist;
        }
//        idB += dist/disc_len;
        idB++;
//        delta_id = (dist-max_dist)/dist;
//        if (delta_id<1) delta_id = 1;
//        idB += delta_id;

    }
    return closest_id;
}


bool MCS_BranchWinder::direction_is_forward(int idA, int idB) {
    if (closed) {
        return (pmod(idB-idA,num_bp) <= pmod(idA-idB,num_bp));
    }
    return (idA<idB);
}




