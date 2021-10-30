#include "MCS_CSrot_2dpot.h"


MCS_CSrot_2dpot::MCS_CSrot_2dpot(Chain * ch,const std::vector<long long int> & seedseq, int selrange_min, int selrange_max)
: MCStep(ch,seedseq)
{
    move_name = "MCS_CSrot_2dpot";
    /*
    Initializes the controling variables of the rotational Crankshaft move
    */

    kz = 5;
    kz = kz/2;

    fac  = 0.25;
//    fac  = 0.1;

    /*
        fac for plasmids
    */
    sigma  = fac*sqrt(chain->get_T()/300)*std::sqrt(disc_len/0.34);
    closed = chain->topology_closed();

    // set the range of the random generators selecting the hinges
    if (selrange_max >= num_bp/2)    selrange_max = num_bp/2;
    if (selrange_min < 2)            selrange_min = 2;
    if (selrange_max < selrange_min) selrange_max = selrange_min;

    if (closed) {
        decltype(genhinge1.param()) new_range(0, num_bp-1);
        genhinge1.param(new_range);

        rge_id1  = 0;
        rge_id2  = num_bp-1;
        rge_span = rge_id2-rge_id1;
    }
    else {
        int rfrom, rto;
        if (!chain->fixed_first_orientation())  rfrom = 0;
        else                                    rfrom = 1;
        rto   = num_bp-1-selrange_min;
        decltype(genhinge1.param()) new_range(rfrom, rto);
        genhinge1.param(new_range);

        rge_id1  = rfrom;
        rge_id2  = num_bp-1;
        rge_span = rge_id2-rge_id1;
    }
    decltype(genhingedist.param()) new_range(selrange_min, selrange_max);
    genhingedist.param(new_range);

////  OLD  //////////////////////////////////////////////////////////////////////
//    if (selrange_max >= num_bp/2)    selrange_max = num_bp/2;
//    if (selrange_min < 2)            selrange_min = 2;
//    if (selrange_max < selrange_min) selrange_max = selrange_min;
//
//    if (closed) {
//        decltype(genhinge1.param()) new_range(0, num_bp-1);
//        genhinge1.param(new_range);
//    }
//    else {
//        int rfrom, rto;
//        if (!chain->fixed_first_orientation())  rfrom = 0;
//        else                                    rfrom = 1;
//        rto   = num_bp-1-selrange_min;
//
//        decltype(genhinge1.param()) new_range(rfrom, rto);
//        genhinge1.param(new_range);
//    }
//    decltype(genhingedist.param()) new_range(selrange_min, selrange_max);
//    genhingedist.param(new_range);
////////////////////////////////////////////////////////////////////////////////

    // intervals of moved segments for exluded volume interactions
    moved_intervals.push_back({-1,-1,1});
    moved_intervals.push_back({0,-1,0});
    moved_intervals.push_back({-1,num_bp-1,0});

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


MCS_CSrot_2dpot::MCS_CSrot_2dpot(Chain * ch,const std::vector<long long int> & seedseq, int selrange_min, int selrange_max, int range_id1, int range_id2)
: MCStep(ch,seedseq)
{
    move_name = "MCS_CSrot_2dpot";
    /*
    Initializes the controling variables of the rotational Crankshaft move
    */

    kz = 5;
    kz = kz/2;

    fac  = 0.25;
//    fac  = 0.1;

    /*
        fac for plasmids
    */
    sigma  = fac*sqrt(chain->get_T()/300)*std::sqrt(disc_len/0.34);
    closed = chain->topology_closed();

    // set the range of the random generators selecting the hinges
    if (selrange_max >= num_bp/2)    selrange_max = num_bp/2;
    if (selrange_min < 2)            selrange_min = 2;
    if (selrange_max < selrange_min) selrange_max = selrange_min;

    rge_id1  = pmod(range_id1  ,num_bp);
    rge_id2  = pmod(range_id2+1,num_bp);

    if (closed) {
        int range = pmod(rge_id2-rge_id1,num_bp);
        if (range == num_bp || range == 0) {
            rge_id1  = 0;
            rge_id2  = num_bp-1;
            decltype(genhinge1.param()) new_range(rge_id1,rge_id2);
            genhinge1.param(new_range);
        }
        else {
            closed_topol_restricted_range = true;
            int upper = larger(rge_id1,rge_id1+range-selrange_min);
            decltype(genhinge1.param()) new_range(rge_id1, upper);
            genhinge1.param(new_range);
        }
    }
    else {
        if (rge_id1 < 0 || rge_id1 > num_bp-2) {
            std::cout << "Error invalid Range selection in MCS_CSrot_2dpot. For open chains range_id1 needs to be within the range [0,num_bp-2]!" << std::endl;
            std::exit(0);
        }
        if (rge_id2 < 0 || rge_id2 > num_bp-1) {
            std::cout << "Error invalid Range selection in MCS_CSrot_2dpot. For open chains range_id2 needs to be within the range [0,num_bp-2]!" << std::endl;
            std::exit(0);
        }
        if (rge_id2 < rge_id1) {
            std::cout << "Error invalid Range selection in MCS_CSrot_2dpot. For open chains range_id2 needs to be larger or equal to range_id1!" << std::endl;
            std::exit(0);
        }
        if (chain->fixed_first_orientation() && rge_id1 == 0) {
            rge_id1 = 1;
        }

        int upper = larger(rge_id1,rge_id2-selrange_min);
        decltype(genhinge1.param()) new_range(rge_id1, upper);
        genhinge1.param(new_range);
    }
    rge_span = pmod(rge_id2-rge_id1,num_bp);

    decltype(genhingedist.param()) new_range(selrange_min, selrange_max);
    genhingedist.param(new_range);

    // intervals of moved segments for exluded volume interactions
    moved_intervals.push_back({-1,-1,1});
    moved_intervals.push_back({0,-1,0});
    moved_intervals.push_back({-1,num_bp-1,0});

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


MCS_CSrot_2dpot::~MCS_CSrot_2dpot() {

}

void MCS_CSrot_2dpot::update_settings() {
    sigma = fac*sqrt(chain->get_T()/300)*std::sqrt(disc_len/0.34);
}


//bool MCS_CSrot_2dpot::MC_move() {
//    int idA = genhinge1(gen);
//    int idB = genhingedist(gen) + idA;
//    return rotation(idA,idB);
//}

bool MCS_CSrot_2dpot::MC_move() {

    int             idA,idB;
    int             idAn,idBn;
    int             hingedist;
    double          theta;
    arma::mat       Rlab,TA_rot,TBn_rot;
	arma::colvec    Theta;
    double          deltaE = 0;

    idA       = genhinge1(gen);
    hingedist = genhingedist(gen);
    idB       = idA + hingedist;
    if (closed) {
        idA    = pmod(idA,num_bp);
        idAn   = pmod(idA-1,num_bp);
        idB    = pmod(idB,num_bp);
        if (closed_topol_restricted_range){
            int span = pmod(idB-rge_id1,num_bp);
            if (span > rge_span || span < hingedist) {
                // Check to keep idBn within the selected range.
                idB = rge_id2;
                hingedist = pmod(idB-idA,num_bp);
            }
        }
        idBn   = pmod(idB-1,num_bp);
    }
    else {
        idAn = idA-1;
        if (idB > rge_id2 ) {
            idB = rge_id2;
            hingedist = idB-idA;
        }
        idBn = idB-1;
    }

    // Trial Move
    theta  = normaldist(gen)*sigma;

    Theta = pos->col(idB) - pos->col(idA);
    Theta = Theta / arma::norm(Theta) * theta;

    Rlab = getRotMat(Theta);

    TA_rot 	= Rlab*triads->slice(idA);
    TBn_rot = Rlab*triads->slice(idBn);

    BPS[idBn]->propose_move(TBn_rot,triads->slice(idB));

    if (idAn >= 0) {
        BPS[idAn]->propose_move(triads->slice(idAn),TA_rot);
        deltaE += BPS[idAn]->eval_delta_energy();

        changed_bps = {idAn,idBn};
    }
    else {
        changed_bps = {idBn};
    }

    deltaE  += BPS[idBn]->eval_delta_energy();

    /*
        Confining Potential
    */
    int mvd = pmod(idB-idA-1,num_bp);
    double z_Eold = 0;
    double z_Enew = 0;
    int k;
    arma::colvec v,vp;
    arma::colvec pnew;

    for (int i=0;i<mvd;i++) {
        k = pmod(idA+1+i,num_bp);
        v = pos->col(k)-pos->col(idA);
        vp = Rlab*v;
        pnew = pos->col(idA) + vp;
        z_Eold += kz*pos->col(k)(2)*pos->col(k)(2);
        z_Enew += kz*pnew(2)*pnew(2);
    }

    deltaE += z_Enew-z_Eold;

    if (exp(-deltaE) <= uniformdist(gen)) {
		return false;
    }
    // If Move is accepted rotate full segments
	else {
        triads->slice(idA)  = TA_rot;
        triads->slice(idBn) = TBn_rot;

        int id,idm1;
        for (int i=idA+1;i<(idA+hingedist-1);i++) {
            id   = pmod( i   ,num_bp);
            idm1 = pmod((i-1),num_bp);

            triads->slice(id) = Rlab*triads->slice(id);
            pos->col(id)      = pos->col(idm1) + triads->slice(idm1).col(2)*disc_len;
        }
        idm1 = pmod(idBn-1,num_bp);
        pos->col(idBn) = pos->col(idm1) + triads->slice(idm1).col(2)*disc_len;

        //////////////////////////////////////////////
        // Assign moved intervals
        int a  = pmod(idA,num_bp);
        int b  = idBn;
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
    }
    return true;
}


bool MCS_CSrot_2dpot::rotation(int idA, int idB) {

//    int             idA,idB;
    int             idAn,idBn;
    double          theta;
    arma::mat       Rlab,TA_rot,TBn_rot;
	arma::colvec    Theta;
    double          deltaE = 0;

//    idA = genhinge1(gen);
//    idB = genhingedist(gen) + idA;
    if (closed) {
        idA    = pmod(idA  ,num_bp);
        idAn   = pmod(idA-1,num_bp);
        idB    = pmod(idB  ,num_bp);
        idBn   = pmod(idB-1,num_bp);
    }
    else {
        idAn = idA-1;
        if (idB >= num_bp ) {
            idB = num_bp-1;
        }
        idBn = idB-1;
    }

    // Trial Move
    theta  = normaldist(gen)*sigma;

    Theta = pos->col(idB) - pos->col(idA);
    Theta = Theta / arma::norm(Theta) * theta;

    Rlab = getRotMat(Theta);

    TA_rot 	= Rlab*triads->slice(idA);
    TBn_rot = Rlab*triads->slice(idBn);

    BPS[idBn]->propose_move(TBn_rot,triads->slice(idB));

    if (idAn >= 0) {
        BPS[idAn]->propose_move(triads->slice(idAn),TA_rot);
        deltaE += BPS[idAn]->eval_delta_energy();

        changed_bps = {idAn,idBn};
    }
    else {
        changed_bps = {idBn};
    }

    deltaE  += BPS[idBn]->eval_delta_energy();

    bool accepted = (exp(-deltaE) > uniformdist(gen));
    // If Move is accepted rotate full segments
	if (accepted) {
        triads->slice(idA)  = TA_rot;
        triads->slice(idBn) = TBn_rot;

        int id,idm1;
        int elems;
        if (idBn < idA) {
            elems = idBn-idA+num_bp;
        }
        else {
            elems = idBn-idA;
        }
        for (int i=idA+1;i<(idA+elems);i++) {
            id   = pmod(i,num_bp);
            idm1 = pmod((i-1),num_bp);

            triads->slice(id) = Rlab*triads->slice(id);
            pos->col(id)      = pos->col(idm1) + triads->slice(idm1).col(2)*disc_len;
        }
        idm1 = pmod(idBn-1,num_bp);
        pos->col(idBn) = pos->col(idm1) + triads->slice(idm1).col(2)*disc_len;
    }

    for (unsigned i=0;i<changed_bps.n_elem;i++) {
        BPS[changed_bps(i)]->set_move(accepted);
    }

    return accepted;
}


