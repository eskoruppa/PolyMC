#include "MCS_Slither2d.h"




MCS_Slither2d::MCS_Slither2d(Chain * ch, double seg_size_min, double seg_size_max, arma::colvec normal_to_plane)
: MCStep(ch),normal(normal_to_plane)
{
    move_name = "MCS_Sli2d";

    /*
    Initializes the controling variables of the rotational Crankshaft move
    */
    fac    = 0.01;
    sigma  = fac*sqrt(chain->get_T()/300)*std::sqrt(disc_len/0.34);

    closed  = chain->topology_closed();
    normal  = normal/arma::norm(normal);
    max_rot = 0.2;

    // set the range of the random generators selecting the hinges
    int n_beads_seg_min = seg_size_min/disc_len;
    int n_beads_seg_max = seg_size_max/disc_len;

    if (n_beads_seg_max >= num_bp/8)  n_beads_seg_max = num_bp/8;
    if (n_beads_seg_min < 1)          n_beads_seg_min = 1;

    cout << "Slider4: Minimum Segment Size = " << n_beads_seg_min << endl;
    cout << "Slider4: Maximum Segment Size = " << n_beads_seg_max << endl;


    if (closed) {
        decltype(genhingeA.param()) new_range(0, num_bp-1);
        genhingeA.param(new_range);
    }
    else {
        int rfrom, rto;
        rfrom = 1;
        rto   = num_bp-1-n_beads_seg_min*3;

        decltype(genhingeA.param()) new_range(rfrom, rto);
        genhingeA.param(new_range);
    }
    decltype(genhingedist.param()) new_range(n_beads_seg_min, n_beads_seg_max);
    genhingedist.param(new_range);

    // intervals of moved segments for exluded volume interactions
    moved_intervals.push_back({-1,-1,-1});
    moved_intervals.push_back({-1,-1,-1});
    moved_intervals.push_back({-1,-1,-1});
    moved_intervals.push_back({-1,-1,-1});
    moved_intervals.push_back({-1,-1,-1});

    changed_bps = arma::zeros<arma::ivec>(4);

    suitablility_key[MCS_FORCE_ACTIVE]              = true;  // force_active
    suitablility_key[MCS_TORQUE_ACTIVE]             = true;  // torque_active
    suitablility_key[MCS_TOPOLOGY_CLOSED]           = true;  // topology_closed
    suitablility_key[MCS_FIXED_LINK]                = true;  // link_fixed
    suitablility_key[MCS_FIXED_TERMINI]             = true;  // fixed_termini
    suitablility_key[MCS_FIXED_TERMINI_RADIAL]      = true;  // fixed_termini_radial
    suitablility_key[MCS_FIXED_FIRST_ORIENTATION]   = true;  // fixed_first_orientation
    suitablility_key[MCS_FIXED_LAST_ORIENTATION]    = true;  // fixed_last_orientation

    requires_EV_check=true;
}


MCS_Slither2d::~MCS_Slither2d() {

}

void MCS_Slither2d::update_settings() {
    sigma = fac*sqrt(chain->get_T()/300)*std::sqrt(disc_len/0.34);
}


//bool MCS_CSrot::MC_move() {
//    int idA = genhinge1(gen);
//    int idB = genhingedist(gen) + idA;
//    return rotation(idA,idB);
//}

bool MCS_Slither2d::MC_move() {

    int          idA,idB,idC,idD;
    int          idAm1,idBm1,idCm1,idDm1;
    int          hingedist;
    arma::colvec posA,posB,posC,posD;
    arma::colvec v1,v2,v3;
    double       len_v1,len_v2,len_v3;
    double       theta1,theta2,theta3;
//    arma::colvec Theta1,Theta2,Theta3;
    arma::colvec rot1,rot3;
    double       dot1,dot3;
    arma::mat    R_Theta1,R_Theta2,R_Theta3;

    arma::colvec posBp,posCp;
    arma::colvec w,nw;
    double       len_w,x,y;
    arma::colvec wt,nwt;
    double       len_wt;

    arma::colvec v1p,v2p,v3p;
    double       delta_E;

    double rge_before,rge_after;

    /*
        Calculate relevant indices
    */

    hingedist = genhingedist(gen);
    idA   = genhingeA(gen);
    idB   = pmod(idA+hingedist,num_bp);
    idC   = pmod(idB+hingedist,num_bp);
    idD   = pmod(idC+hingedist,num_bp);
    idAm1 = pmod(idA-1,num_bp);
    idBm1 = pmod(idB-1,num_bp);
    idCm1 = pmod(idC-1,num_bp);
    idDm1 = pmod(idD-1,num_bp);

    if (!chain->topology_closed()) {
        bool found = false;
        while (idA+3*hingedist >= num_bp) {
            hingedist = genhingedist(gen);
            idA   = genhingeA(gen);
        }
        idB   = pmod(idA+hingedist,num_bp);
        idC   = pmod(idB+hingedist,num_bp);
        idD   = pmod(idC+hingedist,num_bp);
        idAm1 = pmod(idA-1,num_bp);
        idBm1 = pmod(idB-1,num_bp);
        idCm1 = pmod(idC-1,num_bp);
        idDm1 = pmod(idD-1,num_bp);
    }



    /*
        Asign hinge positions and connection vectors
    */
    posA = pos->col(idA);
    posB = pos->col(idB);
    posC = pos->col(idC);
    posD = pos->col(idD);

    v1 = posB-posA;
    v2 = posC-posB;
    v3 = posD-posC;
    len_v1 = arma::norm(v1);
    len_v2 = arma::norm(v2);
    len_v3 = arma::norm(v3);


    arma::colvec cnct = posD-posA;
//    arma::colvec rge;
//    rge = max_angle(v2, cnct, len_v1+len_v3, max_rot);
//    rge_before = rge(1) - rge(0);
//    theta2   = rge(0) + uniformdist(gen) * rge_before;
    theta2 = -max_rot + uniformdist(gen)*2*max_rot;

//    theta1   = normaldist(gen)*sigma;
    R_Theta2 = getRotMat(theta2*normal);

    v2p   = R_Theta2*v2;

//    if (R_Theta2.has_nan()) {
//        cout << "Nan in R_Theta2" << endl;
//        cout << theta2 << endl;
//        cout << rge_before << endl;
//        cout << rge(0) << endl;
//    }



    w = cnct - v2p;
    len_w = arma::norm(w);


//    if (len_w > len_v1+len_v3) {
//        // REMOVE!
//        cout << "impossible move (large)" << endl;
//        cout << len_w << endl;
//        cout << len_v1+len_v3 << endl;
//        cout << len_v1+len_v3-len_w << endl;
//        cout << theta2 << endl;
//        cout << rge.t();
//        cout << arma::norm(cnct) - len_v1+len_v2+len_v3 << endl;
//        cout << arma::norm(v1+v2+v3) << endl;
//        cout << arma::norm(cnct) << endl;
//        return false;
////        std::exit(0);
//    }
//
//    if (len_w < 1e-3) {
//        // REMOVE!
//        cout << "impossible move (small)" << endl;
//        return false;
////        std::exit(0);
//    }

    if (len_w < 1e-3 || len_w > len_v1+len_v3) {
        /*
            The move is directy rejected if the new configurations would lead to overlapping B and D or
            when the segments v2 and v3 cannot be rearranged to bridge the distance between B and D.
        */
        return false;
    }

    nw     = w/len_w;
    x      = (len_v1*len_v1-len_v3*len_v3+len_w*len_w)/(2*len_w);
    nwt    = arma::cross(normal,nw);
    y      = std::sqrt(len_v1*len_v1-x*x);


//    arma::colvec k = (posB-posA) - arma::dot(posB-posA,nw);
//    k = k/arma::norm(k);
//    if (arma::dot(k,nwt)<=0) {
//        posBp = posA+x*nw+y*nwt;
//    }
//    else {
//        posBp = posA+x*nw-y*nwt;
//    }

    if (uniformdist(gen)<0.5) {
        posBp = posA+x*nw+y*nwt;
    }
    else {
        posBp = posA+x*nw-y*nwt;
    }

    posCp = posBp + v2p;
    v1p   = posBp - posA;
    v3p   = posD  - posCp;

//    if (len_v2 > len_w || len_v3 > len_w) {
//        cout << "Possible problem!!!" << endl;
////        std::exit(0);
//    }


    /*
        Calculate the rotation angle around the planar normal necessary to rotate v1 into v1p and
        v2 into v2p respectively.
    */
    rot1   = arma::cross(v1,v1p)/(len_v1*len_v1);
    dot1   = arma::dot(normal,rot1);



//    if (1-std::abs(dot1)<1e-14) {
//        return false;
//    }
//
//    if (1-std::abs(dot1)<1e-10) {
//        cout << "Too Close to 90 deg rotation (v1)" << endl;
//        cout << "dot   = " << 1-std::abs(dot1) << endl;
//        cout << "range = " << rge.t();
//        cout << "theta = " << theta2 << endl;
//        cout << v1.t()/len_v1;
//        cout << v1p.t()/len_v1;
//
//        cout << "############################" << endl;
//        cout << "############################" << endl;
//        cout << "############################" << endl;
//
//
//        int tot = 1e8;
//        int cnt = 0;
//
//        double theta1s = 0;
//        double theta2s = 0;
//
//
//        arma::colvec k = (posB-posA) - arma::dot(posB-posA,nw);
//        k = k/arma::norm(k);
//
//        for (int i=0;i<tot;i++) {
//
//            theta2   = rge(0) + uniformdist(gen) * rge_before;
//            R_Theta2 = getRotMat(theta2*normal);
//
//            v2p   = R_Theta2*v2;
//            w = cnct - v2p;
//            len_w = arma::norm(w);
//
//            nw     = w/len_w;
//            x      = (len_v1*len_v1-len_v3*len_v3+len_w*len_w)/(2*len_w);
//            nwt    = arma::cross(normal,nw);
//            y      = std::sqrt(len_v1*len_v1-x*x);
//
//
//
//
////            if (arma::dot(k,nwt)<=0) {
////                posBp = posA+x*nw+y*nwt;
////            }
////            else {
////                posBp = posA+x*nw-y*nwt;
////            }
//
//            if (uniformdist(gen)<0.5) {
//                posBp = posA+x*nw+y*nwt;
//            }
//            else {
//                posBp = posA+x*nw-y*nwt;
//            }
//
//            posCp = posBp + v2p;
//            v1p   = posBp - posA;
//            v3p   = posD  - posCp;
//
//            rot1   = arma::cross(v1,v1p)/(len_v1*len_v1);
//            dot1   = arma::dot(normal,rot1);
//            if (1-std::abs(dot1)<1e-14) {
////            if (std::abs(arma::dot(v1,v1p)/(len_v1*len_v1))<1e-3) {
//                cnt++;
//            }
//            else {
//                theta1 = std::asin(dot1);
//                if (arma::dot(v1,v1p)<0) {
//                    theta1 = M_PI-theta1;
//                }
//                theta1s += std::abs(theta1);
//                theta2s += std::abs(theta2);
//            }
//        }
//
//        cout << "90deg occured in " << cnt*1./tot*100 << "% of the trials" << endl;
//        cout << "Mean Theta1 = " << theta1s/(tot-cnt)*180/M_PI << endl;
//        cout << "Mean Theta2 = " << theta2s/(tot-cnt)*180/M_PI << endl;
//
//
//
//
//        return false;
//    }

    theta1 = std::asin(dot1);
    if (theta1 != theta1) {
//        cout << "nan theta1 found! " << endl;
//        cout << 1-std::abs(dot1) << endl;
//        cout << dot1 << endl;
//        cout << rot1.t();
//        cout << v1.t();
//        cout << v1p.t();
//        cout << len_v1 << endl;
        return false;
//        std::exit(0);
    }


    if (arma::dot(v1,v1p)<0) {
        theta1 = M_PI-theta1;
    }
    R_Theta1 = getRotMat(theta1*normal);

    rot3   = arma::cross(v3,v3p)/(len_v3*len_v3);
    dot3   = arma::dot(normal,rot3);
    if (1-std::abs(dot3)<1e-14) {
        return false;
    }
    /*
    if (1-std::abs(dot3)<1e-10) {
        cout << "Too Close to 90 deg rotation (v3)" << endl;
        cout << "dot   = " << 1-std::abs(dot3) << endl;
        cout << "range = " << rge.t();
        cout << "theta = " << theta2 << endl;
        cout << v3.t()/len_v3;
        cout << v3p.t()/len_v3;
        return false;
    }
    */
    theta3 = std::asin(dot3);
    if (theta3 != theta3) {
        cout << "nan theta3 found! " << endl;
        cout << 1-std::abs(dot3) << endl;
        std::exit(0);
    }

    if (arma::dot(v3,v3p)<0) {
        theta3 = M_PI-theta3;
    }
    R_Theta3 = getRotMat(theta3*normal);


    arma::mat TA,TBm1,TB,TCm1,TC,TDm1;

    TA   = R_Theta1*triads->slice(idA);
    TBm1 = R_Theta1*triads->slice(idBm1);
    TB   = R_Theta2*triads->slice(idB);
    TCm1 = R_Theta2*triads->slice(idCm1);
    TC   = R_Theta3*triads->slice(idC);
    TDm1 = R_Theta3*triads->slice(idDm1);

    delta_E = 0;



    BPS[idAm1]->propose_move(triads->slice(idAm1),TA);
    BPS[idBm1]->propose_move(TBm1,TB);
    BPS[idCm1]->propose_move(TCm1,TC);
    BPS[idDm1]->propose_move(TDm1,triads->slice(idD));

    delta_E += BPS[idAm1]->eval_delta_energy();
    delta_E += BPS[idBm1]->eval_delta_energy();
    delta_E += BPS[idCm1]->eval_delta_energy();
    delta_E += BPS[idDm1]->eval_delta_energy();

    changed_bps(0) = idAm1;
    changed_bps(1) = idBm1;
    changed_bps(2) = idCm1;
    changed_bps(3) = idDm1;

    /*
    arma::mat Mnorm = arma::zeros(3,3);
    Mnorm(0,0) = 50;
    Mnorm(1,1) = 50;
    Mnorm(2,2) = 100;

    arma::mat Mflex = arma::zeros(3,3);
    Mflex(0,0) = 2;
    Mflex(1,1) = 2;
    Mflex(2,2) = 1;

    arma::mat M;
    arma::mat OmegaO,OmegaN;

    double dEcheck = 0;

    M = Mnorm;
    if (idAm1==num_bp-1) {
        M = Mflex;
    }
    OmegaO = ExtractTheta(triads->slice(idAm1).t()*triads->slice(idA));
    OmegaN = ExtractTheta(triads->slice(idAm1).t()*TA);
    dEcheck -= 1./(2*disc_len)*arma::dot(OmegaO,M*OmegaO);
    dEcheck += 1./(2*disc_len)*arma::dot(OmegaN,M*OmegaN);

    M = Mnorm;
    if (idBm1==num_bp-1) {
        M = Mflex;
    }
    OmegaO = ExtractTheta(triads->slice(idBm1).t()*triads->slice(idB));
    OmegaN = ExtractTheta(TBm1.t()*TB);
    dEcheck -= 1./(2*disc_len)*arma::dot(OmegaO,M*OmegaO);
    dEcheck += 1./(2*disc_len)*arma::dot(OmegaN,M*OmegaN);

    M = Mnorm;
    if (idCm1==num_bp-1) {
        M = Mflex;
    }
    OmegaO = ExtractTheta(triads->slice(idCm1).t()*triads->slice(idC));
    OmegaN = ExtractTheta(TCm1.t()*TC);
    dEcheck -= 1./(2*disc_len)*arma::dot(OmegaO,M*OmegaO);
    dEcheck += 1./(2*disc_len)*arma::dot(OmegaN,M*OmegaN);

    M = Mnorm;
    if (idDm1==num_bp-1) {
        M = Mflex;
    }
    OmegaO = ExtractTheta(triads->slice(idDm1).t()*triads->slice(idD));
    OmegaN = ExtractTheta(TDm1.t()*triads->slice(idD));
    dEcheck -= 1./(2*disc_len)*arma::dot(OmegaO,M*OmegaO);
    dEcheck += 1./(2*disc_len)*arma::dot(OmegaN,M*OmegaN);


//    cout << "### Energy ###" << endl;
//    cout << "diff  = " << delta_E - dEcheck << endl;
//    cout << "dE    = " << delta_E << endl;
//    cout << "check = " << dEcheck << endl;

    if (std::abs(delta_E - dEcheck) > 1e-12) {
        cout << "problem with energies!" << endl;
        std::exit(0);
    }
    */




//    cout << "B -> A" << endl;
//    rge = max_angle(v2p, cnct, len_v1+len_v3, max_rot);
//    rge_after = rge(1) - rge(0);

//    cout << "rnge = " << rge_before << " -> " << rge_after << endl;

//    delta_E = 0;

//    if (exp(-delta_E) * rge_before/rge_after <= uniformdist(gen)) {
    if (exp(-delta_E)  <= uniformdist(gen)) {
//        cout << "move rejected" << endl;
		return false;
    }

    /*
        Perform the rotations and translations of the full segments
    */

    int id,idp1,num;

    /*
        Segment AB
    */
    idp1               = pmod(idA+1,num_bp);
    triads->slice(idA) = TA;
    pos->col(idp1)     = pos->col(idA) + TA.col(2)*disc_len;

    num  = pmod(idB-idA,num_bp);
    for (int a=1;a<num;a++) {
        id   = idp1;
        idp1 = pmod(id+1,num_bp);
        triads->slice(id) = R_Theta1*triads->slice(id);
        pos->col(idp1)    = pos->col(id) + triads->slice(id).col(2)*disc_len;
    }

    /*
        Segment BC
    */
    idp1               = pmod(idB+1,num_bp);
    triads->slice(idB) = TB;
    pos->col(idp1)     = pos->col(idB) + TB.col(2)*disc_len;

    num  = pmod(idC-idB,num_bp);
    for (int a=1;a<num;a++) {
        id   = idp1;
        idp1 = pmod(id+1,num_bp);
        triads->slice(id) = R_Theta2*triads->slice(id);
        pos->col(idp1)    = pos->col(id) + triads->slice(id).col(2)*disc_len;
    }

    /*
        Segment CD
    */
    idp1               = pmod(idC+1,num_bp);
    triads->slice(idC) = TC;
    pos->col(idp1)     = pos->col(idC) + TC.col(2)*disc_len;

    num  = pmod(idD-idC,num_bp);
    for (int a=1;a<num;a++) {
        id   = idp1;
        idp1 = pmod(id+1,num_bp);
        triads->slice(id) = R_Theta3*triads->slice(id);
        pos->col(idp1)    = pos->col(id) + triads->slice(id).col(2)*disc_len;
    }

    /*
        testing for instability
    */
    double dist = arma::norm(posD-pos->col(idD));
    if (dist>1e-8) {
        cout << "BranchWinder Displacement exceeds 1e-8 (" << dist << ")" << endl;
//        EV->revert_to_backup();
        chain->restore_consistency();
//        std::exit(0);
        return false;
    }
    if (pos->has_nan()) {
        cout << "Nan in pos (sli2d)" << endl;
        cout << "theta1 = " << theta1 << endl;
        cout << "rot1   = " << rot1.t() << endl;
        cout << "dot1   = " << arma::dot(normal,rot1) << endl;
        cout << arma::dot(v1,v1p) << endl;
        cout << R_Theta1;
        cout << "theta2 = " << theta2 << endl;
        cout << R_Theta2;
        cout << "theta3 = " << theta3 << endl;
        cout << "rot3   = " << rot3.t() << endl;
        cout << "dot3   = " << arma::dot(normal,rot3) << endl;
        cout << arma::dot(v3,v3p) << endl;
        cout << R_Theta3;

        cout << "--------------------" << endl;
//        cout << rge.t();
        cout << theta2 << endl;
        std::exit(0);
    }
    if (triads->has_nan()) {
        cout << "Nan in triads (sli2d)" << endl;
        std::exit(0);
    }

    /*
        Assign Moved Intervals
    */

    if (idA < idD) {
        moved_intervals[0](0) = 0;
        moved_intervals[0](1) = idA;
        moved_intervals[0](2) = 0;

        moved_intervals[1](0) = idA+1;
        moved_intervals[1](1) = idB;
        moved_intervals[1](2) = 1000;

        moved_intervals[2](0) = idB+1;
        moved_intervals[2](1) = idC;
        moved_intervals[2](2) = 1001;

        moved_intervals[3](0) = idC+1;
        moved_intervals[3](1) = idD-1;
        moved_intervals[3](2) = 1002;

        moved_intervals[4](0) = idD;
        moved_intervals[4](1) = num_bp-1;
        moved_intervals[4](2) = 0;
    }
    else {
        if (idB < idD) {
            moved_intervals[0](0) = 0;
            moved_intervals[0](1) = idB;
            moved_intervals[0](2) = 1000;

            moved_intervals[1](0) = idB+1;
            moved_intervals[1](1) = idC;
            moved_intervals[1](2) = 1001;

            moved_intervals[2](0) = idC+1;
            moved_intervals[2](1) = idD-1;
            moved_intervals[2](2) = 1002;

            moved_intervals[3](0) = idD;
            moved_intervals[3](1) = idA;
            moved_intervals[3](2) = 0;

            moved_intervals[4](0) = idA;
            moved_intervals[4](1) = num_bp-1;
            moved_intervals[4](2) = 1000;
        }
        else {
            if (idC < idD) {
                moved_intervals[0](0) = 0;
                moved_intervals[0](1) = idC;
                moved_intervals[0](2) = 1001;

                moved_intervals[1](0) = idC+1;
                moved_intervals[1](1) = idD-1;
                moved_intervals[1](2) = 1002;

                moved_intervals[2](0) = idD;
                moved_intervals[2](1) = idA;
                moved_intervals[2](2) = 0;

                moved_intervals[3](0) = idA+1;
                moved_intervals[3](1) = idB;
                moved_intervals[3](2) = 1000;

                moved_intervals[4](0) = idB+1;
                moved_intervals[4](1) = num_bp-1;
                moved_intervals[4](2) = 1001;
            }
            else {
                moved_intervals[0](0) = 0;
                moved_intervals[0](1) = idD-1;
                moved_intervals[0](2) = 1002;

                moved_intervals[1](0) = idD;
                moved_intervals[1](1) = idA;
                moved_intervals[1](2) = 0;

                moved_intervals[2](0) = idA+1;
                moved_intervals[2](1) = idB;
                moved_intervals[2](2) = 1000;

                moved_intervals[3](0) = idB+1;
                moved_intervals[3](1) = idC;
                moved_intervals[3](2) = 1001;

                moved_intervals[4](0) = idC+1;
                moved_intervals[4](1) = num_bp-1;
                moved_intervals[4](2) = 1002;
            }
        }
    }
    return true;
}


arma::colvec MCS_Slither2d::max_angle(arma::colvec& v1, arma::colvec& w, double lv23, double max_select) {
    double lv1,lw;
    lv1 = arma::norm(v1);
    lw  = arma::norm(w);
    arma::colvec nv1 = v1 / lv1;
    arma::colvec nw  = w  / lw;

    double theta_max;
    if (lv23 > lw+lv1) {
        theta_max = M_PI;
    }
    else {
        double max_forwardsq = lw*lw+lv1*lv1;
        if (lv23*lv23 <= max_forwardsq) {
            /*
                Max angle is smaller than or equal to pi/2
            */
            double lw1 = (lv1*lv1-lv23*lv23+lw*lw)/(2*lw);
            theta_max = std::acos(lw1/lv1);
        }
        else {
            /*
                Max angle is larger than pi/2
            */
            double lw1 = (lv23*lv23-lv1*lv1-lw*lw)/(2*lw);
            theta_max = M_PI-std::acos(lw1/lv1);
        }
    }

    double angle_from,angle_to;
    double dot = arma::dot(nv1,nw);
    if (dot > 1) {
        dot = 1;
    }
    if (dot < -1) {
        dot = -1;
    }
    double current_theta = std::acos(dot);
    if (current_theta > theta_max) {
        cout << "Error in max theta calculation" << endl;
        cout << "current_theta = " << current_theta << endl;
        cout << "theta_max     = " << theta_max << endl;
        return {0,0};
    }

    arma::colvec o = arma::cross(normal,nw);
    if (arma::dot(o,nv1)>=0) {
        angle_from = theta_max+current_theta;
        angle_to   = theta_max-current_theta;

    }
    else {
        angle_from = theta_max-current_theta;
        angle_to   = theta_max+current_theta;
    }
    if (angle_from > max_select) {
        angle_from = max_select;
    }
    if (angle_to > max_select) {
        angle_to = max_select;
    }
    if ( angle_from != angle_from ) {
        cout << "angle_from is nan!" << endl;
        cout << theta_max << endl;
        cout << current_theta << endl;
        cout << 1-arma::dot(nv1,nw) << endl;
        cout << arma::dot(v1,w) << endl;
        cout << nv1.t();
        cout << nw.t();
        cout << max_select << endl;
    }
    if ( angle_to != angle_to ) {
        cout << "angle_to is nan!" << endl;
        cout << theta_max << endl;
        cout << current_theta << endl;
        cout << arma::dot(nv1,nw) << endl;
        cout << arma::dot(v1,w) << endl;
        cout << max_select << endl;
    }


    return {-angle_from+1e-8,angle_to-1e-8};
}




