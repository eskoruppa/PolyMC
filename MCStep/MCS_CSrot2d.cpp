#include "MCS_CSrot2d.h"


MCS_CSrot2d::MCS_CSrot2d(Chain * ch,const std::vector<long long int> & seedseq, int selrange_min, int selrange_max,arma::colvec normal)
: MCStep(ch,seedseq),normal(normal)
{
    move_name = "MCS_CSrot2d";
    /*
    Initializes the controling variables of the rotational Crankshaft move
    */
    closed = chain->topology_closed();
    normal = normal/arma::norm(normal);

    // set the range of the random generators selecting the hinges
    selrange_max = selrange_max/2;
    selrange_min = selrange_min/2;

    if (selrange_max >= num_bp/4)    selrange_max = num_bp/4;
    if (selrange_min < 2)            selrange_min = 2;
    if (selrange_max < selrange_min) selrange_max = selrange_min;

    if (closed) {
        decltype(genhinge1.param()) new_range(0, num_bp-1);
        genhinge1.param(new_range);
    }
    else {
        int rfrom, rto;
        if (!chain->fixed_first_orientation())  rfrom = 0;
        else                                    rfrom = 1;
        rto   = num_bp-1-selrange_min;

        decltype(genhinge1.param()) new_range(rfrom, rto);
        genhinge1.param(new_range);
    }
    decltype(genhingedist.param()) new_range(selrange_min, selrange_max);
    genhingedist.param(new_range);

    // intervals of moved segments for exluded volume interactions
    moved_intervals.push_back({-1,-1,1});
    moved_intervals.push_back({0,-1,0});
    moved_intervals.push_back({-1,num_bp-1,0});

    changed_bps = arma::zeros<arma::ivec>(3);

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


MCS_CSrot2d::~MCS_CSrot2d() {

}

void MCS_CSrot2d::update_settings() {
}

bool MCS_CSrot2d::MC_move() {

    int         idA,idB,idC;
    int         idAm1,idBm1,idCm1;
    int         hingedist;
    arma::colvec posA,posB,posC;
    arma::colvec posBp;

    arma::colvec w,nw;
    double       len_w;
    arma::colvec nv;
    double       vAB,wAB;

    arma::colvec v1,v2,v1p,v2p;
    double       len_v1,len_v2;
    arma::colvec rotA,rotB;
    double       thetaA,thetaB;
    arma::mat    R_thetaA,R_thetaB;

    arma::mat   TA,TBm1,TB,TCm1;
    double      delta_E = 0;

    /*
        choose hinge ids
    */
    idA       = genhinge1(gen);
    hingedist = genhingedist(gen);
    idB       = pmod(hingedist+idA,num_bp);
    idC       = pmod(hingedist+idB,num_bp);
    idAm1     = pmod(idA-1,num_bp);
    idBm1     = pmod(idB-1,num_bp);
    idCm1     = pmod(idC-1,num_bp);

    if (!chain->topology_closed()) {
        while (idA+2*hingedist >= num_bp || idA == 0) {
            hingedist = genhingedist(gen);
            idA       = genhinge1(gen);
        }
        idB       = pmod(hingedist+idA,num_bp);
        idC       = pmod(hingedist+idB,num_bp);
        idAm1     = pmod(idA-1,num_bp);
        idBm1     = pmod(idB-1,num_bp);
        idCm1     = pmod(idC-1,num_bp);
//        cout << idA << " " << idB << " " << idC << " " << endl;
    }

    /*
        assign initial hinge positions
    */
    posA = pos->col(idA);
    posB = pos->col(idB);
    posC = pos->col(idC);

    /*
        vector connecting A and C
    */
    w       = posC-posA;
    len_w   = arma::norm(w);
    nw      = w/len_w;

    /*
        second vector in the plane orthogonal to w
    */
    nv = arma::cross(normal,nw);
    nv = nv/arma::norm(nv);

    /*
        Decompose the vector AB (v1) into w and v components. For the new configuration
        the component of v is flipped.
    */
    v1 = posB-posA;
    vAB = arma::dot(v1,nv);
    wAB = arma::dot(v1,nw);

    /*
        Calculate the new position of B and find the old and new vectors connecting
        AB and BC.
    */
    posBp = posA+wAB*nw-vAB*nv;
    v1p   = posBp-posA;
    v2    = posC-posB;
    v2p   = posC-posBp;

    len_v1 = arma::norm(v1);
    len_v2 = arma::norm(v2);

    /*
        Calculate the rotation angle around the planar normal necessary to rotate v1 into v1p and
        v2 into v2p respectively.
    */
    rotA   = arma::cross(v1,v1p)/(len_v1*len_v1);
    thetaA = std::asin(arma::dot(normal,rotA));
    if (arma::dot(v1,v1p)<0) {
        thetaA = M_PI-thetaA;
    }
    R_thetaA = getRotMat(thetaA*normal);

    rotB   = arma::cross(v2,v2p)/(len_v2*len_v2);
    thetaB = std::asin(arma::dot(normal,rotB));
    if (arma::dot(v2,v2p)<0) {
        thetaB = M_PI-thetaB;
    }
    R_thetaB = getRotMat(thetaB*normal);

    if (R_thetaA.has_nan() || R_thetaB.has_nan()) {
        std::cout << "Nan-value in rotation matrix" << std::endl;
        return false;
    }
    if (pos->has_nan()) {
        std::cout << "Nan in pos (rot2d)" << std::endl;
        std::exit(0);
    }
    if (triads->has_nan()) {
        std::cout << "Nan in triads (rot2d)" << std::endl;
        std::exit(0);
    }


    /*
        Calculate the new triads that lead to a change in elastic energy.
    */
    TA   = R_thetaA*triads->slice(idA);
    TBm1 = R_thetaA*triads->slice(idBm1);
    TB   = R_thetaB*triads->slice(idB);
    TCm1 = R_thetaB*triads->slice(idCm1);

    /*
        Calculate change in elastic energy
    */
    BPS[idAm1]->propose_move(triads->slice(idAm1),TA);
    BPS[idBm1]->propose_move(TBm1,TB);
    BPS[idCm1]->propose_move(TCm1,triads->slice(idC));

    delta_E += BPS[idAm1]->eval_delta_energy();
    delta_E += BPS[idBm1]->eval_delta_energy();
    delta_E += BPS[idCm1]->eval_delta_energy();

    changed_bps(0) = idAm1;
    changed_bps(1) = idBm1;
    changed_bps(2) = idCm1;

    /*
        Metropolis Step
    */
    if (exp(-delta_E) <= uniformdist(gen)) {
		return false;
    }

    /*
        displace all positions between idA and idC and rotate the respective triads
    */

    int id,idm1,num;

    /*
        Rotate triads A
    */
    triads->slice(idA) = TA;
    /*
        Segments between A and B
    */
    num = pmod(idB-idA,num_bp);
    id  = idA;
    for (int i=1;i<num;i++) {
        idm1 = id;
        id   = pmod(idm1+1,num_bp);

        pos->col(id)      = pos->col(idm1)+triads->slice(idm1).col(2)*disc_len;
        triads->slice(id) = R_thetaA*triads->slice(id);
    }

    /*
        Segment B
    */
//    if (id != pmod(idB-1,num_bp)) {
//        cout << "wrong ID!" << endl;
//        cout << id << " != " << pmod(idB-1,num_bp) << endl;
//        cout << "idA = " << idA << endl;
//        cout << "idB = " << idB << endl;
//        cout << "num = " << num << endl;
//        std::exit(0);
//    }

    pos->col(idB)      = pos->col(id)+triads->slice(id).col(2)*disc_len;
    triads->slice(idB) = R_thetaB*triads->slice(idB);

    /*
        Check if this positions of B conincides with the prior calculated one.
    */
//    if (arma::norm(pos->col(idB)-posBp)>1e-10) {
//        cout << "problem with posBp! (" << arma::norm(pos->col(idB)-posBp) << ")" << endl;
//        cout << thetaA << endl;
//        cout << num << endl;
//    }
//    if (arma::norm(v1p-R_thetaA*v1) > 1e-10) {
//        cout << "wrong rotation v1!" << endl;
//        std::exit(0);
//    }

    /*
        Segments between B and C
    */
    num = pmod(idC-idB,num_bp);
    id  = idB;
    for (int i=1;i<num;i++) {
        idm1 = id;
        id   = pmod(idm1+1,num_bp);

        pos->col(id)      = pos->col(idm1)+triads->slice(idm1).col(2)*disc_len;
        triads->slice(id) = R_thetaB*triads->slice(id);
    }

    /*
        Check check if the chain integrity is preserved
    */
//    arma::colvec newC = pos->col(idCm1)+triads->slice(idCm1).col(2)*disc_len;
//    if (arma::norm(pos->col(idC)-newC)>1e-10) {
//        cout << "problem with posC! (" << arma::norm(pos->col(idC)-newC) << ")" << endl;
//        cout << thetaB << endl;
//        cout << num << endl;
//        cout << idA << " " << idC << endl;
//        std::exit(0);
//    }
//    if (arma::norm(v2p-R_thetaB*v2) > 1e-10) {
//        cout << "wrong rotation v2!" << endl;
//        std::exit(0);
//    }

    //////////////////////////////////////////////
    // Assign moved intervals
    int a  = pmod(idA,num_bp);
    int b  = idCm1;
    if (a<b) {
        if (a==0) {
            if (b==num_bp-1) {
                // situation V
//                    cout << "SITUATION V" << endl;
                // moved interval spans full range
                moved_intervals[0](EV_FROM) = a;
                moved_intervals[0](EV_TO)   = b;
                moved_intervals[0](EV_TYPE) = 1000;

                moved_intervals[1](EV_TYPE) = -1;
                moved_intervals[2](EV_TYPE) = -1;
            }
            else {
                // situation II
//                    cout << "SITUATION II" << endl;
                // moved interval starting at first bp
                moved_intervals[0](EV_FROM) = a;
                moved_intervals[0](EV_TO)   = b;
                moved_intervals[0](EV_TYPE) = 1000;

                moved_intervals[1](EV_FROM) = b+1;
                moved_intervals[1](EV_TO)   = num_bp-1;
                moved_intervals[1](EV_TYPE) = 0;

                moved_intervals[2](EV_TYPE) = -1;
            }
        }
        else {
            if (b==num_bp-1) {
                // situation III
//                    cout << "SITUATION III" << endl;
                // moved interval ending at last bp
                moved_intervals[1](EV_FROM) = a;
                moved_intervals[1](EV_TO)   = b;
                moved_intervals[1](EV_TYPE) = 1000;

                moved_intervals[0](EV_FROM) = 0;
                moved_intervals[0](EV_TO)   = a-1;
                moved_intervals[0](EV_TYPE) = 0;

                moved_intervals[2](EV_TYPE) = -1;
            }
            else {
                // situation I
//                    cout << "SITUATION I" << endl;
                // moved interval in the middle
                moved_intervals[1](EV_FROM) = a;
                moved_intervals[1](EV_TO)   = b;
                moved_intervals[1](EV_TYPE) = 1000;

                moved_intervals[0](EV_FROM) = 0;
                moved_intervals[0](EV_TO)   = a-1;
                moved_intervals[0](EV_TYPE) = 0;

                moved_intervals[2](EV_FROM) = b+1;
                moved_intervals[2](EV_TO)   = num_bp-1;
                moved_intervals[2](EV_TYPE) = 0;
            }
        }
    }
    else {
        // situation IV
//            cout << "SITUATION IV" << endl;
//            cout << idA << " " << idB << endl;
//            cout << a << " " << b << endl;
        // moved interval crosses boundary
        moved_intervals[0](EV_FROM) = 0;
        moved_intervals[0](EV_TO)   = b;
        moved_intervals[0](EV_TYPE) = 1000;

        moved_intervals[2](EV_FROM) = a;
        moved_intervals[2](EV_TO)   = num_bp-1;
        moved_intervals[2](EV_TYPE) = 1000;

        moved_intervals[1](EV_FROM) = b+1;
        moved_intervals[1](EV_TO)   = a-1;
        moved_intervals[1](EV_TYPE) = 0;

    }
    return true;
}

