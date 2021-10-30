#include "MCS_Slither2d.h"

/*
    UNFINSHED!
    Potentially doesn't satisfy detailed balance.
*/


MCS_Slither2d::MCS_Slither2d(Chain * ch,const std::vector<long long int> & seedseq, double seg_size_min, double seg_size_max, arma::colvec normal_to_plane, bool rotate_tail)
: MCStep(ch,seedseq),normal(normal_to_plane), rot_tail(rotate_tail)
{
    move_name = "MCS_Sli2d";

    /*
    Initializes the controling variables of the rotational Crankshaft move
    */
    fac    = 0.01;
    sigma  = fac*sqrt(chain->get_T()/300)*std::sqrt(disc_len/0.34);

    closed  = chain->topology_closed();
    normal  = normal/arma::norm(normal);
    max_rot = 0.05;

    // set the range of the random generators selecting the hinges
    int n_beads_seg_min = seg_size_min/disc_len;
    int n_beads_seg_max = seg_size_max/disc_len;

    if (n_beads_seg_max >= num_bp/8)  n_beads_seg_max = num_bp/8;
    if (n_beads_seg_min < 1)          n_beads_seg_min = 1;

    std::cout << "Slider4: Minimum Segment Size = " << n_beads_seg_min << std::endl;
    std::cout << "Slider4: Maximum Segment Size = " << n_beads_seg_max << std::endl;


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
    if (!rot_tail) moved_intervals.push_back({-1,-1,-1});

    int num_changed;
    if (rot_tail)   num_changed=3;
    else            num_changed=4;
    changed_bps = arma::zeros<arma::ivec>(num_changed);

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

bool MCS_Slither2d::MC_move() {

    int          idA,idB,idC,idD;
    int          idAm1,idBm1,idCm1,idDm1;
    int          hingedist;
    arma::colvec posA,posB,posC,posD;
    arma::colvec v1,v2,v3;
//    double       len_v1;
    double       len_v2,len_v3;
    double       theta1,theta2,theta3;
    arma::colvec rot1,rot2,rot3;
    arma::mat    R_Theta1,R_Theta2,R_Theta3;

    arma::colvec posBp,posCp;
    arma::colvec w,nw;
    double       len_w,len_w2;
    arma::colvec wt,nwt;
    double       len_wt;

    arma::colvec v1p,v2p,v3p;
    double       delta_E;

    double rge_before;

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
    len_v2 = arma::norm(v2);
    len_v3 = arma::norm(v3);

    arma::colvec cnct = posD-posA;
    arma::colvec rge;
//    rge = max_angle(v1, cnct, len_v2+len_v3, max_rot);
    rge = max_angle(v1, cnct, len_v2+len_v3, 2*M_PI);
    rge_before = rge(1) - rge(0);


//    theta1   = rge(0) + uniformdist(gen) * rge_before;
    theta1 = (uniformdist(gen)-0.5)*0.01;


    R_Theta1 = getRotMat(theta1*normal);

    v1p   = R_Theta1*v1;
    posBp = posA+v1p;

    w     = posD - posBp;
    len_w = arma::norm(w);

    if (len_w < 1e-3 || len_w > len_v2+len_v3) {
        /*
            The move is directy rejected if the new configurations would lead to overlapping B and D or
            when the segments v2 and v3 cannot be rearranged to bridge the distance between B and D.
        */
        return false;
    }

    nw     = w/len_w;
    len_w2 = (len_v2*len_v2-len_v3*len_v3+len_w*len_w)/(2*len_w);
    nwt    = arma::cross(normal,nw);
    len_wt = std::sqrt(len_v2*len_v2-len_w2*len_w2);



    if (uniformdist(gen)<0.5) {
        posCp = posBp+len_w2*nw+len_wt*nwt;
    }
    else {
        posCp = posBp+len_w2*nw-len_wt*nwt;
    }

    v2p = posCp - posBp;
    v3p = posD  - posCp;

    /*
        Calculate the rotation angle around the planar normal necessary to rotate v1 into v1p and
        v2 into v2p respectively.
    */
    rot2   = arma::cross(v2,v2p)/(len_v2*len_v2);
    theta2 = std::asin(arma::dot(normal,rot2));
    if (theta2 != theta2) {
        return false;
    }

    if (arma::dot(v2,v2p)<0) {
        theta2 = M_PI-theta2;
    }
    R_Theta2 = getRotMat(theta2*normal);

    rot3   = arma::cross(v3,v3p)/(len_v3*len_v3);
    theta3 = std::asin(arma::dot(normal,rot3));
    if (theta3 != theta3) {
        return false;
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
    if (!rot_tail) BPS[idDm1]->propose_move(TDm1,triads->slice(idD));

    delta_E += BPS[idAm1]->eval_delta_energy();
    delta_E += BPS[idBm1]->eval_delta_energy();
    delta_E += BPS[idCm1]->eval_delta_energy();
    if (!rot_tail) delta_E += BPS[idDm1]->eval_delta_energy();

    changed_bps(0) = idAm1;
    changed_bps(1) = idBm1;
    changed_bps(2) = idCm1;
    if (!rot_tail) changed_bps(3) = idDm1;


    arma::colvec nx,ny;
    nx = posD-posA;
    nx = nx/arma::norm(nx);
    ny = arma::cross(normal,nx);
    double fac = cal_gfrac(v1,v2,v3,v1p,v2p,v3p,nx,ny,rge_before);

//    if (exp(-delta_E) <= uniformdist(gen)) {
    if (exp(-delta_E) * fac <= uniformdist(gen)) {
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
        Rotate Tail
    */
    if (rot_tail) {
        triads->slice(idD) = R_Theta3*triads->slice(idD);
        for (int i=idD+1;i<num_bp;i++) {
            pos->col(i)    = pos->col(i-1) + triads->slice(i-1).col(2)*disc_len;
            triads->slice(i) = R_Theta3*triads->slice(i);
        }
    }

    /*
        testing for instability
    */
    double dist = arma::norm(posD-pos->col(idD));
    if (dist>1e-8) {
        std::cout << "Slither2d Displacement exceeds 1e-8 (" << dist << ")" << std::endl;
//        EV->revert_to_backup();
        std::exit(0);
        return false;
    }
    if (pos->has_nan()) {
        std::cout << "Nan in pos (sli2d)" << std::endl;
        std::cout << "theta1 = " << theta1 << std::endl;
        std::cout << R_Theta1;
        std::cout << "theta2 = " << theta2 << std::endl;
        std::cout << "rot2   = " << rot2.t() << std::endl;
        std::cout << "dot2   = " << arma::dot(normal,rot2) << std::endl;
        std::cout << arma::dot(v2,v2p) << std::endl;
        std::cout << R_Theta2;
        std::cout << "theta3 = " << theta3 << std::endl;
        std::cout << "rot3   = " << rot3.t() << std::endl;
        std::cout << "dot3   = " << arma::dot(normal,rot3) << std::endl;
        std::cout << arma::dot(v3,v3p) << std::endl;
        std::cout << R_Theta3;

        std::cout << "--------------------" << std::endl;
//        std::cout << rge.t();
        std::cout << theta2 << std::endl;
        std::exit(0);
    }
    if (triads->has_nan()) {
        std::cout << "Nan in triads (sli2d)" << std::endl;
        std::exit(0);
    }

    /*
        Assign Moved Intervals
    */
    if (!rot_tail) {
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
    }
    else {
        /*
            If the tail is rotated
        */
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
        moved_intervals[3](1) = num_bp-1;
        moved_intervals[3](2) = 1002;
    }
    return true;
}


arma::colvec MCS_Slither2d::max_angle(arma::colvec& v1, arma::colvec& w, double lv23, double max_select) {
    double lv1,lw;
    lv1 = arma::norm(v1);
    lw  = arma::norm(w);
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

    arma::colvec nv1 = v1 / lv1;
    arma::colvec nw  = w  / lw;
    double angle_from,angle_to;
    double current_theta = std::acos(arma::dot(nv1,nw));
    if (current_theta > theta_max) {
        std::cout << "Error in max theta calculation" << std::endl;
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
    return {-angle_from,angle_to};
}

double MCS_Slither2d::cal_dens(arma::colvec& u, arma::colvec& v, arma::colvec& w, arma::colvec& nx, arma::colvec& ny, double rho1) {
    double u1,u2,v1,w1,w2; //,v2
    u1 = arma::dot(u,nx);
    u2 = arma::dot(u,ny);
    v1 = arma::dot(v,nx);
//    v2 = arma::dot(v,ny);
    w1 = arma::dot(w,nx);
    w2 = arma::dot(w,ny);

    double denom = (v1*w2+w1*w2+w1*u2);
    double dT2 = (w1*u2-w2*u1)/denom;
    double dT4 = (u1*w2+u1*u2+v1*u2)/denom;
    double dT3 = -dT4;

    double dens = 1./std::sqrt(1+dT2*dT2+dT3*dT3);
    return dens;
}

double MCS_Slither2d::cal_gfrac(arma::colvec& u, arma::colvec& v, arma::colvec& w,arma::colvec& up, arma::colvec& vp, arma::colvec& wp, arma::colvec& nx, arma::colvec& ny, double rge) {
    double rho1 = std::abs(1./rge);

    double rhoBA = cal_dens(u,v,w,nx,ny,rho1);
    double rhoAB = cal_dens(up,vp,wp,nx,ny,rho1);

    return rhoBA/rhoAB;
}




