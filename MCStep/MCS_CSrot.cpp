#include "MCS_CSrot.h"

#define MCS_CSROT_FAC 0.25


MCS_CSrot::MCS_CSrot(Chain * ch,const std::vector<long long int> & seedseq, int selrange_min, int selrange_max)
: MCStep(ch,seedseq)
{
    move_name = "MCS_CSrot";
    /*
    Initializes the controling variables of the rotational Crankshaft move
    */
    fac  = MCS_CSROT_FAC;
    /*
        fac for plasmids
    */
    sigma  = fac*sqrt(chain->get_T()/300)*std::sqrt(disc_len/0.34);
    closed = chain->topology_closed();

    // set the range of the random generators selecting the hinges
    if (closed) {
        if (selrange_max >= num_bp/2)    selrange_max = num_bp/2;
    }
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

    requires_EV_check   = true;
    requires_pair_check = true;
}


MCS_CSrot::MCS_CSrot(Chain * ch,const std::vector<long long int> & seedseq, int selrange_min, int selrange_max, int range_id1, int range_id2)
: MCStep(ch,seedseq)
/*

    If the chain is open range_id2 cannot be larger than num_bp-2
*/
{
    move_name = "MCS_CSrot";
    /*
    Initializes the controling variables of the rotational Crankshaft move
    */
    fac  = MCS_CSROT_FAC;
    /*
        fac for plasmids
    */
    sigma  = fac*sqrt(chain->get_T()/300)*std::sqrt(disc_len/0.34);
    closed = chain->topology_closed();

    // set the range of the random generators selecting the hinges
    if (closed) {
        if (selrange_max >= num_bp/2)    selrange_max = num_bp/2;
    }
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
            std::cout << "Error invalid Range selection in MCS_CSrot. For open chains range_id1 needs to be within the range [0,num_bp-2]!" << std::endl;
            std::exit(0);
        }
        if (rge_id2 < 0 || rge_id2 > num_bp-1 || range_id2+1 == num_bp) {
            std::cout << "Error invalid Range selection in MCS_CSrot. For open chains range_id2 needs to be within the range [0,num_bp-2]!" << std::endl;
            std::exit(0);
        }
        if (rge_id2 < rge_id1) {
            std::cout << "Error invalid Range selection in MCS_CSrot. For open chains range_id2 needs to be larger or equal to range_id1!" << std::endl;
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

    requires_EV_check   = true;
    requires_pair_check = true;
}



MCS_CSrot::~MCS_CSrot() {

}

void MCS_CSrot::update_settings() {
    sigma = fac*sqrt(chain->get_T()/300)*std::sqrt(disc_len/0.34);
}


//bool MCS_CSrot::MC_move() {
//    int idA = genhinge1(gen);
//    int idB = genhingedist(gen) + idA;
//    return rotation(idA,idB);
//}

bool MCS_CSrot::MC_move() {

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
    theta = normaldist(gen)*sigma;
    Theta = pos->col(idB) - pos->col(idA);
    Theta = Theta / arma::norm(Theta) * theta;

    Rlab    = getRotMat(Theta);
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
        int a  = idA; //pmod(idA+1,num_bp);
        int b  = idBn;
//        std::cout << " CSrot:    " << a << " - " << b << std::endl;
        if (a<b) {
            if (a==0) {
                if (b==num_bp-1) {
                    // situation V
                    // moved interval spans full range
                    moved_intervals[0](EV_FROM) = a;
                    moved_intervals[0](EV_TO)   = b;
                    moved_intervals[0](EV_TYPE) = 1000;

                    moved_intervals[1](EV_TYPE) = -1;
                    moved_intervals[2](EV_TYPE) = -1;
                }
                else {
                    // situation II
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


bool MCS_CSrot::rotation(int idA, int idB) {

//    int             idA,idB;gen_trial_conf
    int             idAn,idBn;
    double          theta;
    arma::mat       Rlab,TA_rot,TBn_rot;
	arma::colvec    Theta;
    double          deltaE = 0;

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
    theta = normaldist(gen)*sigma;
    Theta = pos->col(idB) - pos->col(idA);
    Theta = Theta / arma::norm(Theta) * theta;

    Rlab    = getRotMat(Theta);
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


double MCS_CSrot::gen_trial_conf() {

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
    theta = normaldist(gen)*sigma;
    Theta = pos->col(idB) - pos->col(idA);
    Theta = Theta / arma::norm(Theta) * theta;

    Rlab    = getRotMat(Theta);
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

//    if (exp(-deltaE) <= uniformdist(gen)) {
//		return false;
//    }
//    // If Move is accepted rotate full segments
//	else {

	/*
        Modification of the configuration
	*/

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
    int a  = idA; //pmod(idA+1,num_bp);
    int b  = idBn;
//        std::cout << " CSrot:    " << a << " - " << b << std::endl;
    if (a<b) {
        if (a==0) {
            if (b==num_bp-1) {
                // situation V
                // moved interval spans full range
                moved_intervals[0](EV_FROM) = a;
                moved_intervals[0](EV_TO)   = b;
                moved_intervals[0](EV_TYPE) = 1000;

                moved_intervals[1](EV_TYPE) = -1;
                moved_intervals[2](EV_TYPE) = -1;
            }
            else {
                // situation II
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
    return deltaE;
}



