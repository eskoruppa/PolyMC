#include "MCS_Slither4.h"


MCS_Slither4::MCS_Slither4(Chain * ch,const std::vector<long long int> & seedseq, double seg_size_min, double seg_size_max)
: MCStep(ch,seedseq)
{
    move_name = "MCS_Slith";

    /*
    Initializes the controling variables of the rotational Crankshaft move
    */

    double fac_temp = std::sqrt(chain->get_T()/300);

    fac_theta = 0.18;
    fac_phi   = 0.02;

    fac_theta = 0.15;
    fac_phi   = 0.018;

    arma::mat covmat = chain->get_avg_cov();
    sigma_theta    = arma::zeros(3);
    sigma_theta(0)  = fac_theta*sqrt(covmat(0,0));
    sigma_theta(1)  = fac_theta*sqrt(covmat(1,1));
    sigma_theta(2)  = fac_theta*sqrt(covmat(2,2));


    sigma_phi = fac_phi*fac_temp;

    closed = chain->topology_closed();

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


MCS_Slither4::~MCS_Slither4() {

}

void MCS_Slither4::update_settings() {
    double fac_temp = std::sqrt(chain->get_T()/300);
    arma::mat covmat = chain->get_avg_cov();
    sigma_theta(0)  = fac_theta*sqrt(covmat(0,0));
    sigma_theta(1)  = fac_theta*sqrt(covmat(1,1));
    sigma_theta(2)  = fac_theta*sqrt(covmat(2,2));
    sigma_phi = fac_phi*fac_temp;
}


//bool MCS_CSrot::MC_move() {
//    int idA = genhinge1(gen);
//    int idB = genhingedist(gen) + idA;
//    return rotation(idA,idB);
//}

bool MCS_Slither4::MC_move() {

    int idA,idB,idC,idD;
    int idAm1,idBm1,idCm1,idDm1;
    int hingedist;
    arma::colvec posA,posB,posC,posD;
    arma::colvec posBp,posCp;
    arma::colvec v1,v2,v3;
    double len_v2,len_v3,len_w;

    arma::colvec Theta;
    arma::mat    R_Theta;

    bool forward_move;

    arma::colvec v1p,v2p,v3p;
    arma::colvec w,nw;
    double len_w2,len_w3;

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

//    if (uniformdist(gen)<0.5) {
    if (true) {
        forward_move = true;

        posA = pos->col(idA);
        posB = pos->col(idB);
        posC = pos->col(idC);
        posD = pos->col(idD);

        v1 = posB-posA;
        v2 = posC-posB;
        v3 = posD-posC;

        Theta    = {normaldist(gen)*sigma_theta(0),normaldist(gen)*sigma_theta(1),normaldist(gen)*sigma_theta(2)};
        R_Theta  = getRotMat(triads->slice(idAm1)*Theta);
    }
    else {
        forward_move = false;

        posA = pos->col(idD);
        posB = pos->col(idC);
        posC = pos->col(idB);
        posD = pos->col(idA);

        v1 = posB-posA;
        v2 = posC-posB;
        v3 = posD-posC;

        Theta    = {normaldist(gen)*sigma_theta(0),normaldist(gen)*sigma_theta(1),normaldist(gen)*sigma_theta(2)};
        R_Theta  = getRotMat(triads->slice(idD)*Theta);
    }

    v1p   = R_Theta*v1;
    posBp = posA+v1p;

    w      = posD-posBp;
    len_w  = arma::norm(w);
    len_v2 = arma::norm(v2);
    len_v3 = arma::norm(v3);


    if ((len_v2+len_v3-len_w)<1e-12) {
        if (len_w > len_v2+len_v3) {
            /*
                If the segments v2 and v3 cannot bridge the distance between posBp and posD
                the move is directly rejected. There cannot be any retries for finding an
                acceptable rotation because that would violate the Metropolis coniditons of
                symmetric selection probabilites.
            */
            return false;
        }
        nw     = w/len_w;
        posCp  = nw*len_v2;
    }
    else {
        arma::colvec r2,r3,nr2,nr3;
        double len_r2,len_r3;
        arma::colvec posM;

        nw     = w/len_w;
        len_w2 = (len_w*len_w+len_v2*len_v2-len_v3*len_v3)/(2*len_w);
        len_w3 = len_w-len_w2;

        posM = posBp + len_w2*nw;

        r2     = v2-arma::dot(nw,v2)*nw;
        len_r2 = arma::norm(r2);
        nr2    = r2/len_r2;

        r3     = -(v3-arma::dot(nw,v3)*nw);
        len_r3 = arma::norm(r3);
        nr3    = r3/len_r3;

        double phi23, delta_phi;
        arma::colvec rotax;
        double circ_radius;
        double dotr2r3;

        dotr2r3 = arma::dot(nr2,nr3);

        if (std::abs(dotr2r3) < 1e-8) {
            phi23 = 0;
            rotax = nw;
        }
        else {
            rotax = arma::cross(nr2,nr3);
            phi23 = std::asin(arma::norm(rotax));
            if (dotr2r3<0) {
                phi23 = M_PI-phi23;
            }
        }

        delta_phi   = normaldist(gen)*sigma_phi;
        circ_radius = std::sqrt(len_v2*len_v2-len_w2*len_w2);
        posCp = posM + getRotMat((0.5*phi23 + delta_phi)*rotax/arma::norm(rotax))*nr2*circ_radius;

    }

    /*

    */
    arma::mat Rv1,Rv2,Rv3;
    double len_v2p,len_v3p;
    arma::colvec nv2,nv3,nv2p,nv3p;

    v2p = posCp-posBp;
    len_v2p = arma::norm(v2p);
    v3p = posD -posCp;
    len_v3p = arma::norm(v3p);

    if ( std::abs(len_v2 - len_v2p)>1e-11 || std::abs(len_v3 - len_v3p)>1e-11  ) {
        std::cout << "BranchWinder: Potential Instability: std::abs(len_v2 - len_v2p)>1e-11 || std::abs(len_v3 - len_v3p)>1e-11" << std::endl;
        return false;
    }

    nv2p = v2p/len_v2p;
    nv3p = v3p/len_v3p;
    nv2 = v2/len_v2;
    nv3 = v3/len_v3;

    arma::colvec axis;
    double len_axis;
    double psi;

    axis = arma::cross(nv2,nv2p);
    len_axis = arma::norm(axis);
    psi = std::asin(len_axis);
    if (arma::dot(nv2,nv2p)<0) {
        psi = M_PI-psi;
    }
    Rv2 = getRotMat(psi*axis/len_axis);


    axis = arma::cross(nv3,nv3p);
    len_axis = arma::norm(axis);
    psi = std::asin(len_axis);
    if (arma::dot(nv3,nv3p)<0) {
        psi = M_PI-psi;
    }
    Rv3 = getRotMat(psi*axis/len_axis);


    arma::mat RA,RB,RC;


    /*
        Put the rotation matrices in ascending order along the actual chain
    */
    if (forward_move) {
        RA = R_Theta;
        RB = Rv2;
        RC = Rv3;
    }
    else {
        RA = Rv3;
        RB = Rv2;
        RC = R_Theta;
    }


    /*
        Evaluate change in energy
    */

    arma::mat TA,TBm1,TB,TCm1,TC,TDm1;

    TA   = RA*triads->slice(idA);
    TBm1 = RA*triads->slice(idBm1);
    TB   = RB*triads->slice(idB);
    TCm1 = RB*triads->slice(idCm1);
    TC   = RC*triads->slice(idC);
    TDm1 = RC*triads->slice(idDm1);

    double delta_E = 0;

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

    if (exp(-delta_E) <= uniformdist(gen)) {
		return false;
    }


    /*
        Perform the rotations and translations of the full segments
    */

    int id,idp1;
    int num;

    /*
        Segment AB
    */
    idp1 = pmod(idA+1,num_bp);

    triads->slice(idA) = TA;
    pos->col(idp1) = pos->col(idA) + TA.col(2)*disc_len;

    num  = pmod(idB-idA,num_bp);
    for (int a=1;a<num;a++) {
        id = idp1;
        idp1 = pmod(id+1,num_bp);
        triads->slice(id) = RA*triads->slice(id);
        pos->col(idp1)    = pos->col(id) + triads->slice(id).col(2)*disc_len;
    }

    /*
        Segment BC
    */
    idp1 = pmod(idB+1,num_bp);

    triads->slice(idB) = TB;
    pos->col(idp1) = pos->col(idB) + TB.col(2)*disc_len;

    num  = pmod(idC-idB,num_bp);
    for (int a=1;a<num;a++) {
        id = idp1;
        idp1 = pmod(id+1,num_bp);
        triads->slice(id) = RB*triads->slice(id);
        pos->col(idp1)    = pos->col(id) + triads->slice(id).col(2)*disc_len;
    }

    /*
        Segment CD
    */
    idp1 = pmod(idC+1,num_bp);

    triads->slice(idC) = TC;
    pos->col(idp1) = pos->col(idC) + TC.col(2)*disc_len;

    num  = pmod(idD-idC,num_bp);
    for (int a=1;a<num;a++) {
        id = idp1;
        idp1 = pmod(id+1,num_bp);
        triads->slice(id) = RC*triads->slice(id);
        pos->col(idp1)    = pos->col(id) + triads->slice(id).col(2)*disc_len;
    }

    /*
        testing for instability
    */
    double dist = arma::norm(posD-pos->col(idD));
    if (dist>1e-10) {
        std::cout << "BranchWinder Displacement exceeds 1e-10 (" << dist << ")" << std::endl;
        EV->revert_to_backup();
        return false;
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






















