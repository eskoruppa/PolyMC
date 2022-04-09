#include "ExcludedVolume.h"

ExVol::ExVol(Chain * ch, double EV_distance) :
chain(ch),
BPS(*ch->get_BPS()),
bp_pos(ch->get_bp_pos()),
triads(ch->get_triads()),
num_bp(ch->get_num_bp()),
num_bps(ch->get_num_bps()),
disc_len(ch->get_disc_len()),
closed_topology(chain->topology_closed()),
counter(0),
counter_reject(0),
counter_RP_reject(0),
repulsion_plane_active(false),
EV_dist(EV_distance),
bp_pos_backup(*chain->get_bp_pos()),
triads_backup(*chain->get_triads())
{


    if (EV_dist < disc_len) {
        throw std::invalid_argument("ExVol::ExVol(): Radius of excluded volume needs to be larger than discretization length!");
    }

    ///////////////////////////////////////////
    // Init EV_beads
    num_bp_per_EV    = int((double)EV_dist/disc_len);
    upper_shift      = num_bp_per_EV-1;
    eff_size_EV_bead = num_bp_per_EV*disc_len;
    num_EV           = int(ceil((double)num_bp/num_bp_per_EV));
    EV_beads         = arma::zeros<arma::ivec>(num_EV);
    for (int i=0;i<num_EV;i++) {
        EV_beads(i) = i*num_bp_per_EV;
    }
    /*
        If there is sufficient overlap between the excluded regions of neighboring EV_beads,
        i.e. if next nearest neighbors overlap for a 90deg angle between the two intermediate
        tangents, excluded volume checks are skipped for the first two neighbours rather than
        just the nearest neighbours.
    */
    neighbour_skip = 1;
    if (2*disc_len*disc_len*num_bp_per_EV*num_bp_per_EV < EV_dist*EV_dist) {
        neighbour_skip = 2;
    }
    neighbour_skip_plus_one          = neighbour_skip+1;
    neighbour_skip_boundary          = neighbour_skip;
    neighbour_skip_boundary_plus_one = neighbour_skip_boundary+1;

    /*
        If the total number of bp is not a multiple of the number of bp per EV_bead and the
        chain has a closed topology, i.e. the first bead is connected to the last, overlap
        between beads displaced by less than (neighbour_skip+1), which can occure around the
        periodic boundary, will be discarded. This is done to prevent false overlap detections.
        Instead, additional static pair checks will be conducted to prevent potential crossings
        in this region.
    */

//    std::cout << neighbour_skip << std::endl;
//    std::cout << num_bp%num_bp_per_EV << std::endl;

    additional_boundcheck=false;
    if (num_bp%num_bp_per_EV!=0 && chain->topology_closed()) {
        std::cout << "mismatch" << std::endl;
        additional_boundcheck=true;
        for (int i=0;i<=neighbour_skip;i++) {
            addboundpairs.push_back({i*num_bp_per_EV,num_bp-(neighbour_skip+1-i)*num_bp_per_EV});
            std::cout << addboundpairs[i].t();
        }
        neighbour_skip_boundary++;
        neighbour_skip_boundary_plus_one++;
    }

    std::cout << std::endl;
    std::cout << "######################################" << std::endl;
    std::cout << "#### INITIALIZING EXCLUDED VOLUME ####" << std::endl;
    std::cout << "######################################" << std::endl;
    std::cout << " Excluded Volume Beads: " << std::endl;
    std::cout << "   number of EV beads: " << num_EV << std::endl;
    std::cout << "   bp per EV bead:     " << num_bp_per_EV << std::endl;
    std::cout << "   Effective size:     " << eff_size_EV_bead << std::endl;
    std::cout << "   Exclusion distance: " << EV_dist << std::endl;
    std::cout << "######################################" << std::endl;
    std::cout << std::endl;
}

ExVol::~ExVol() {

}

////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////// GETTERS AND SETTERS  /////////////////////////////////////////////////////


double ExVol::get_EV_dist() {
    return EV_dist;
}

int    ExVol::get_bp_per_EV() {
    return num_bp_per_EV;
}

arma::mat*  ExVol::get_bp_pos_backup() {
    return &bp_pos_backup;
}

arma::cube* ExVol::get_triads_backup() {
    return &triads_backup;
}

void ExVol::set_backup_conf(arma::mat * backup_pos, arma::cube * backup_triads) {
    bp_pos_backup = *backup_pos;
    triads_backup = *backup_triads;
}

////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////

//inline int ExVol::bp2EV(int bp_id) {
///*
//    returns the id of the EV bead that contains the bp with this bp_id
//*/
//    return bp_id/num_bp_per_EV;
//}

inline int ExVol::bp2EV_upper(int bp_id) {
/*
    returns the id of the first EV bead that has a bp index higher or
    equal bp_id
*/
    return (bp_id+upper_shift)/num_bp_per_EV;
}

inline int ExVol::bp2EV_lower(int bp_id) {
/*
    returns the id of the first EV bead that has a bp index smaller or
    equal bp_id
*/
    return bp_id/num_bp_per_EV;
}

arma::ivec ExVol::interval_bp2EV(const arma::ivec& bp_interval) {
/*
    Transforms a bp interval into an EV interval. This interval contains all the EV
    beads for which the reference bp (by which the positions of the EV is defined)
    lies inside the given bp interval.
*/
    return {bp2EV_upper(bp_interval(0)),bp2EV_lower(bp_interval(1)),bp_interval(2)};
}

inline bool ExVol::within_EV_typeA(int id) {
    return (id == EV_TYPE_A);
}
inline bool ExVol::within_EV_typeB(int id) {
    return (id == EV_TYPE_B);
}
inline bool ExVol::within_EV_typeC(int id) {
    return (EV_TYPE_C_FROM <= id && id <= EV_TYPE_C_TO);
}
inline bool ExVol::within_EV_typeD(int id) {
    return (EV_TYPE_D_FROM <= id && id <= EV_TYPE_D_TO);
}
inline bool ExVol::within_EV_typeE(int id) {
    return (id == EV_TYPE_E);
}


void ExVol::transfer_config(const std::vector<arma::ivec>* moved, arma::mat* bp_pos_from, arma::mat* bp_pos_to, arma::cube* triads_from, arma::cube* triads_to ) {
/*
    Copying the full arma::mat and arma::cube is much faster as counted per element as copying partial mat/cube.
    Hence if the amount of moved elements exceeds 0.2*num_bp the full config will be copied
*/
    *bp_pos_to = *bp_pos_from;
    *triads_to = *triads_from;

//    int from,to;
//    int amount_moved=0;
//    for (int i=0;i<moved->size();i++) {
//        if ((*moved)[i](2) > 0) {
//            from = (*moved)[i](0);
//            to   = (*moved)[i](1);
//            amount_moved+=to-from+1;
//        }
//    }
//
//    if (amount_moved > num_bp*0.2) {
//        *bp_pos_to = *bp_pos_from;
//        *triads_to = *triads_from;
//    }
//    else {
//        for (int i=0;i<moved->size();i++) {
//            if ((*moved)[i](2) > 0) {
//                from = (*moved)[i](0);
//                to   = (*moved)[i](1);
//                if (from <= to) {
//                    bp_pos_to->cols(from,to) = bp_pos_from->cols(from,to);
//                    triads_to->slices(from,to) = triads_from->slices(from,to);
//                }
//            }
//        }
//    }
}

void ExVol::revert_to_backup() {
/*
    Reverts the current chain positions and triads back to the stored backup.
*/
    *bp_pos = bp_pos_backup;
    *triads = triads_backup;
}

void ExVol::set_current_as_backup() {
/*
    Sets the current chain positions and triads as backup.
*/
    bp_pos_backup = *bp_pos;
    triads_backup = *triads;
}

double ExVol::rejection_rate() {
    if (counter>0) {
        return (double)counter_reject/(double)counter;
    }
    else {
        return 0;
    }
}

bool ExVol::check(const std::vector<arma::ivec>* moved) {
    /*
        TODO: REPULSION PLANES
        - use ProximityClusters
        - integrate plane proximity to ProximityClusters
    */
    counter++;

    #ifdef DEBUG_EXVOL
    check_moved_interval_order(moved);
    #endif

    std::vector<arma::ivec> EV_typeA;
    std::vector<arma::ivec> EV_typeB;
    std::vector<arma::ivec> EV_typeC;
    std::vector<arma::ivec> EV_typeD;
    std::vector<arma::ivec> EV_typeE;
    cal_EV_intervals(moved,EV_typeA,EV_typeB,EV_typeC,EV_typeD,EV_typeE);

    bool check = check_intervals(EV_typeA,EV_typeB,EV_typeC,EV_typeD,EV_typeE);
    if (repulsion_plane_active && check) {
//        check = RP_interval(0,num_EV-2);
//        This function doesn't work for now
//        check = RP_check_intervals(EV_typeB,EV_typeC,EV_typeD,EV_typeE);

        if (first_last == REPULSIONPLANE_FIRST_AND_LAST) {
            check = RP_interval(0,num_EV-2);
        }
        else if (first_last == REPULSIONPLANE_ONLY_FIRST) {
            check = RP_interval_first(0,num_EV-2);
        }
        else {
            check = RP_interval_last(0,num_EV-2);
        }
    }

    #ifdef DEBUG_EXVOL
    if (check) {
        if (check_overlap()) {
            std::cout << "Bead Overlap!" << std::endl;
            std::cout << "Moved Intervals:" << std::endl;
            std::cout << "Type A" << std::endl;
            for (int i=0;i<EV_typeA.size();i++) {
                std::cout << (EV_typeA)[i].t();
            }
            std::cout << "Type B" << std::endl;
            for (int i=0;i<EV_typeB.size();i++) {
                std::cout << (EV_typeB)[i].t();
            }
            std::cout << "Type C" << std::endl;
            for (int i=0;i<EV_typeC.size();i++) {
                std::cout << (EV_typeC)[i].t();
            }
            std::cout << "Type D" << std::endl;
            for (int i=0;i<EV_typeD.size();i++) {
                std::cout << (EV_typeD)[i].t();
            }
            std::cout << "Type E" << std::endl;
            for (int i=0;i<EV_typeE.size();i++) {
                std::cout << (EV_typeE)[i].t();
            }
            std::exit(1);
        }
    }
    #endif

    if (check) {
//        std::cout << "accepted!" << std::endl;
        /*
            EV respected. Copy config to backup
        */
        transfer_config(moved, bp_pos, &bp_pos_backup, triads, &triads_backup);
    }
    else {
//        std::cout << "violated!" << std::endl;
        counter_reject++;
        /*
            EV violated. Copy backup to config
        */
        transfer_config(moved, &bp_pos_backup, bp_pos, &triads_backup,triads);
    }

    /*
        Every EV_FULLCHECK_VIOLATION steps carefully check for bead overlaps to
        prevent beads from getting trapped.
    */
    #ifdef EV_FULLCHECK_VIOLATION
    if (counter%EV_FULLCHECK_VIOLATION == 0) {
        if (check_overlap()) {
            /*
                TODO: Add appropriate action
            */
            std::cout << "Simulation stopped due to Excluded volume violation!" << std::endl;
            std::exit(1);
        }
    }
    #endif
    return check;
}



void ExVol::cal_EV_intervals(   const std::vector<arma::ivec>* moved,
                                std::vector<arma::ivec>& EV_typeA,
                                std::vector<arma::ivec>& EV_typeB,
                                std::vector<arma::ivec>& EV_typeC,
                                std::vector<arma::ivec>& EV_typeD,
                                std::vector<arma::ivec>& EV_typeE) {
/*
    Converts the bp intervals into EV intervals and corrects for potential overlap
    - It is critical that the intervals in 'moved' are ordered from left to right.
    TODO: remove the requirement for the moved intervals to be ordered.
*/

    // conversion
    std::vector<arma::ivec> EV_intervals;
    for (unsigned i=0;i<moved->size();i++) {
        if ((*moved)[i](EV_TYPE) >= 0 && (*moved)[i](EV_FROM) <= (*moved)[i](EV_TO)) {
            arma::ivec new_EV_interval = interval_bp2EV((*moved)[i]);
            if (new_EV_interval(0) <= new_EV_interval(1) ) {
                EV_intervals.push_back(new_EV_interval);
            }
        }
    }

    // sort intervals
    for (unsigned i=0;i<EV_intervals.size();i++) {
        if (EV_intervals[i](EV_FROM) <= EV_intervals[i](EV_TO)) {
            if ( within_EV_typeA( EV_intervals[i](EV_TYPE) )) {
                EV_typeA.push_back(EV_intervals[i]);
                continue;
            }
            if ( within_EV_typeB( EV_intervals[i](EV_TYPE) )) {
                EV_typeB.push_back(EV_intervals[i]);
                continue;
            }
            if ( within_EV_typeC( EV_intervals[i](EV_TYPE) )) {
                EV_typeC.push_back(EV_intervals[i]);
                continue;
            }
            if ( within_EV_typeD( EV_intervals[i](EV_TYPE) )) {
                EV_typeD.push_back(EV_intervals[i]);
                continue;
            }
            if ( within_EV_typeE( EV_intervals[i](EV_TYPE) )) {
                EV_typeE.push_back(EV_intervals[i]);
                continue;
            }
        }
    }

    #ifdef DEBUG_EXVOL
    // test consistency
    for (int i=1;i<EV_intervals.size();i++) {
        if (EV_intervals[i](EV_FROM) != EV_intervals[i-1](EV_TO)+1) {
            std::cout << "intervals inconsistent!"  << std::endl;
            for (int j=0;j<EV_intervals.size();j++) {
                std::cout << EV_intervals[j].t();
            }
            std::cout << std::endl;
            for (int i=0;i<moved->size();i++) {
                std::cout << (*moved)[i].t();
            }
            break;
        }
    }
    if (EV_intervals[EV_intervals.size()-1](1) != num_EV-1) {
        std::cout << "intervals inconsistent! (last)" << std::endl;
        std::cout << "Index of last EV bead: " << num_EV-1 << std::endl;
        for (int j=0;j<EV_intervals.size();j++) {
            std::cout << EV_intervals[j].t();
        }
        std::cout << std::endl;
        for (int i=0;i<moved->size();i++) {
            std::cout << (*moved)[i].t();
        }
    }
    #endif
}


bool ExVol::check_intervals(const std::vector<arma::ivec>& EV_typeA,
                            const std::vector<arma::ivec>& EV_typeB,
                            const std::vector<arma::ivec>& EV_typeC,
                            const std::vector<arma::ivec>& EV_typeD,
                            const std::vector<arma::ivec>& EV_typeE) {

    // Check EV_typeC
    for (unsigned tC=0;tC<EV_typeC.size();tC++) {
        // Check with all EV_typeA with singleMove check
        for (int tA=EV_typeA.size()-1;tA>=0;tA--) {
//        for (unsigned tA=0;tA<EV_typeA.size();tA++) {
            if (!check_interval_singleMove(EV_typeC[tC](EV_FROM),EV_typeC[tC](EV_TO),EV_typeA[tA](EV_FROM),EV_typeA[tA](EV_TO))) {
                return false;
            }
        }
        // Check with all EV_typeB with doubleMove check
        for (unsigned tB=0;tB<EV_typeB.size();tB++) {
            if (!check_interval_doubleMove(EV_typeC[tC](EV_FROM),EV_typeC[tC](EV_TO),EV_typeB[tB](EV_FROM),EV_typeB[tB](EV_TO))) {
                return false;
            }
        }
        // Check with all EV_typeD with doubleMove check
        for (unsigned tD=0;tD<EV_typeD.size();tD++) {
            if (!check_interval_doubleMove(EV_typeC[tC](EV_FROM),EV_typeC[tC](EV_TO),EV_typeD[tD](EV_FROM),EV_typeD[tD](EV_TO))) {
                return false;
            }
        }
        // Check with all EV_typeE with doubleMove check
        for (unsigned tE=0;tE<EV_typeE.size();tE++) {
            if (!check_interval_doubleMove(EV_typeC[tC](EV_FROM),EV_typeC[tC](EV_TO),EV_typeE[tE](EV_FROM),EV_typeE[tE](EV_TO))) {
                return false;
            }
        }
        // Check with all other EV_typeC
        for (unsigned tC2=tC+1;tC2<EV_typeC.size();tC2++) {
            /*
                A single interval can be split into two intervals when it crosses the periodic boundary. For typeC intervals
                these intervals should not be mutually checked
            */
            if ( EV_typeC[tC](EV_TYPE) != EV_typeC[tC2](EV_TYPE) ) {
                if (!check_interval_doubleMove(EV_typeC[tC](EV_FROM),EV_typeC[tC](EV_TO),EV_typeC[tC2](EV_FROM),EV_typeC[tC2](EV_TO))) {
                    return false;
                }
            }
        }
    }

    // Check EV_typeD
    for (unsigned tD=0;tD<EV_typeD.size();tD++) {
        // Check with all EV_typeA with singleMove check
        for (unsigned tA=0;tA<EV_typeA.size();tA++) {
            if (!check_interval_singleMove(EV_typeD[tD](EV_FROM),EV_typeD[tD](EV_TO),EV_typeA[tA](EV_FROM),EV_typeA[tA](EV_TO))) {
                return false;
            }
        }
        // Check with all EV_typeB with doubleMove check
        for (unsigned tB=0;tB<EV_typeB.size();tB++) {
            if (!check_interval_doubleMove(EV_typeD[tD](EV_FROM),EV_typeD[tD](EV_TO),EV_typeB[tB](EV_FROM),EV_typeB[tB](EV_TO))) {
                return false;
            }
        }
        // Check with all EV_typeE with doubleMove check
        for (unsigned tE=0;tE<EV_typeE.size();tE++) {
            if (!check_interval_doubleMove(EV_typeD[tD](EV_FROM),EV_typeD[tD](EV_TO),EV_typeE[tE](EV_FROM),EV_typeE[tE](EV_TO))) {
                return false;
            }
        }

        // Check with all other EV_typeD
        for (unsigned tD2=tD+1;tD2<EV_typeD.size();tD2++) {
            /*
                A single interval can be split into two intervals when it crosses the periodic boundary. For typeD intervals
                these intervals SHOULD be mutally checked regardless!
            */
            if (!check_interval_doubleMove(EV_typeD[tD](EV_FROM),EV_typeD[tD](EV_TO),EV_typeD[tD2](EV_FROM),EV_typeD[tD2](EV_TO))) {
                return false;
            }
        }

        // Check within the interval itself
        if (!check_within_interval(EV_typeD[tD](EV_FROM),EV_typeD[tD](EV_TO))) {
            return false;
        }
    }

    // Check EV_typeB
    for (unsigned tB=0;tB<EV_typeB.size();tB++) {
        // Check with all EV_typeA with singleMove check
        for (unsigned tA=0;tA<EV_typeA.size();tA++) {
            if (!check_interval_singleMove(EV_typeB[tB](EV_FROM),EV_typeB[tB](EV_TO),EV_typeA[tA](EV_FROM),EV_typeA[tA](EV_TO))) {
                return false;
            }
        }
        // Check with all EV_typeE with doubleMove check
        for (unsigned tE=0;tE<EV_typeE.size();tE++) {
            if (!check_interval_doubleMove(EV_typeB[tB](EV_FROM),EV_typeB[tB](EV_TO),EV_typeE[tE](EV_FROM),EV_typeE[tE](EV_TO))) {
                return false;
            }
        }
    }

    /*
        Check additional Boundary pairs. These will always be checked regardless on whether the constituent monomers
        are within one of the moved intervals.
    */
    if (additional_boundcheck) {
        double dist;
        for (unsigned i=0;i<addboundpairs.size();i++) {
            dist = doubleMove(addboundpairs[i](0),addboundpairs[i](1));
            if (dist < EV_dist) {
                return false;
            }
        }
    }

    return true;
}

///////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////
////// CHECK INTERVAL SINGLE MOVE


#ifndef EXVOL_USE_OLD_CHECK_INTERVAL
bool ExVol::check_interval_singleMove(int A1, int A2, int B1, int B2) {
/*
    Checks for projected overlap between EV beads in the two intervals limited by
    A1,A2 and B1,B2 respectively. This method assumes that the beads in the B
    interval have not moved, such that only the linear displacement of the beads
    in the A interval have to be considered. The checking is done by the
    singleMove method.
*/

    #ifdef DEBUG_EXVOL
    if (!(A2<B1 || B2 < A1)) {
        std::cout << "Check Intervals " << A1 << " " << A2 << " " << B1 << " " << B2 << std::endl;
        throw std::invalid_argument("ExVol::check_interval_singleMove(): Invalid Intervals! The invervals are not allowed to cross the periodic boundary!");
    }
    #endif

    double dist;
    int a1,a2;

    if (closed_topology) {
    /*
    ///////////////////////////////////////////////////////////////////////////////////////////////
        PERIODIC BOUNDARY CONDITION
    */

        /*
            Left to right order: A - B
        */
        if (A2<B1) {
            a1 = larger(A1,B2-num_EV+neighbour_skip_boundary_plus_one);
            a2 = smaller(A2,B1-neighbour_skip_plus_one);
            /*
                Check for overlaps within the boundary region of A and region B. Pairs within neighbour skip region are omitted.
            */
            int a = A1;
            while (a<=A2) {
                if ( a1<=a && a<=a2) {
                    a=a2+1;
                    continue;
                }
                int b  = larger (B1,a+neighbour_skip_plus_one);
                int b2 = smaller(B2,a+num_EV-neighbour_skip_boundary_plus_one);
                while (b<=b2) {
                    dist = singleMove(EV_beads(a),EV_beads(b));
                    if (dist < EV_dist) {
                        return false;
                    }
                    //b++;
                    b += (dist-EV_dist)/eff_size_EV_bead+1;
                }
                a++;
            }
        }
        /*
            Left to right order: B - A
        */
        else {
            a1 = larger(A1,B2+neighbour_skip_plus_one);
            a2 = smaller(A2,B1+num_EV-neighbour_skip_boundary_plus_one);
            /*
                Check for overlaps within the boundary region of A and region B. Pairs within neighbour skip region are omitted.
            */
            int a = A1;
            while (a<=A2) {
                if ( a1<=a && a<=a2) {
                    a=a2+1;
                    continue;
                }
                int b  = larger (B1,a-num_EV+neighbour_skip_boundary_plus_one);
                int b2 = smaller(B2,a-neighbour_skip_plus_one);
                while (b<=b2) {
                    dist = singleMove(EV_beads(a),EV_beads(b));
                    if (dist < EV_dist) {
                        return false;
                    }
                    //b++;
                    b += (dist-EV_dist)/eff_size_EV_bead+1;
                }
                a++;
            }
        }
    }
    else {
    /*
    ///////////////////////////////////////////////////////////////////////////////////////////////
        NON-PERIODIC BOUNDARY CONDITION
    */
        /*
            Left to right order: A - B
        */
        if (A2<B1) {
            /*
                Define the boundary regions in which pairchecking potentially has to be omitted due to proximity along the chain.
                In this case this region is only at the right boundary of the interval A.
            */
            a1 = A1;
            a2 = smaller(A2,B1-neighbour_skip_plus_one);
            /*
                Check for overlaps within this boundary region of A and region B. Pairs within neighbour skip region are omitted
            */
            for (int a=a2+1;a<=A2;a++) {
                int b = larger(B1,a+neighbour_skip_plus_one);
                while (b<=B2) {
                    dist = singleMove(EV_beads(a),EV_beads(b));
                    if (dist < EV_dist) {
                        return false;
                    }
                    //b++;
                    b += (dist-EV_dist)/eff_size_EV_bead+1;
                }
            }

        }
        /*
            Left to right order: B - A
        */
        else {
            /*
                Define the boundary regions in which pairchecking potentially has to be omitted due to proximity along the chain.
                In this case this region is only at the left boundary of the interval A.
            */
            a1 = larger(A1,B2+neighbour_skip_plus_one);
            a2 = A2;
            /*
                Check for overlaps within this boundary region of A and region B. Pairs within neighbour skip region are omitted
            */
            for (int a=A1;a<a1;a++) {
                int b  = B1;
                int b2 = smaller(B2,a-neighbour_skip_plus_one);
                while (b<=b2) {
                    dist = singleMove(EV_beads(a),EV_beads(b));
                    if (dist < EV_dist) {
                        return false;
                    }
                    //b++;
                    b += (dist-EV_dist)/eff_size_EV_bead+1;
                }
            }

        }
    }
    /*
    ///////////////////////////////////////////////////////////////////////////////////////////////
        Finally check all pairs outside the range of potential neighbour skips. I.e. the bulk of the intervals.
    */
    for (int a=a1;a<=a2;a++) {
        int b = B1;
        while (b<=B2) {
            dist = singleMove(EV_beads(a),EV_beads(b));
            if (dist < EV_dist) {
                return false;
            }
            //b++;
            b += (dist-EV_dist)/eff_size_EV_bead+1;
        }
    }
    return true;
}
#endif


#ifdef EXVOL_USE_OLD_CHECK_INTERVAL
bool ExVol::check_interval_singleMove(int A1, int A2, int B1, int B2) {
/*
    Checks for projected overlap between EV beads in the two intervals limited by
    A1,A2 and B1,B2 respectively. This method assumes that the beads in the B
    interval have not moved, such that only the linear displacement of the beads
    in the A interval have to be considered. The checking is done by the
    singleMove method.
*/

    #ifdef DEBUG_EXVOL
    if (!(A2<B1 || B2 < A1)) {
        std::cout << "Check Intervals " << A1 << " " << A2 << " " << B1 << " " << B2 << std::endl;
        assert(A1<=A2 && B1<=B2);
    }
    if (!(A2<B1 || B2 < A1)) {
        std::cout << "Check Intervals " << A1 << " " << A2 << " " << B1 << " " << B2 << std::endl;
        assert(A2<B1 || B2 < A1);
    }
    #endif

    double dist;
    int b,b1,b2,a1,a2;
    int diff;

    /*
    All pairs outside the range of potential neighbour skips.
    */

    a1 = A1+neighbour_skip;
    if (A1==0 && additional_boundcheck) {
        a1++;
    }
    a2 = A2-neighbour_skip;
    if (A2==num_EV-1 && additional_boundcheck) {
        a2--;
    }

    for (int a=a1;a<=a2;a++) {
        b = B1;
        while (b<=B2) {
            dist = singleMove(EV_beads(a),EV_beads(b));
            if (dist < EV_dist) {
                return false;
            }
//            b++;
            b += (dist-EV_dist)/eff_size_EV_bead+1;
        }
    }

    /*
    The first elements in the interval for which there are potential neighbour skips.
    */
    for (int a=A1;a<a1;a++) {
        int ns=neighbour_skip;
        if (B2>A1) {
            b2 = B2-num_EV;
            if (additional_boundcheck) ns++;
        }
        else {
            b2 = B2;
        }
        int diff = a-b2;
        if (diff <= ns) {
            b2 = a-ns-1;
        }
        b2 = pmod(b2,num_EV);
        b = B1;
        while (b<=b2) {
            dist = singleMove(EV_beads(a),EV_beads(b));
            if (dist < EV_dist) {
                return false;
            }
//            b++;
            b += (dist-EV_dist)/eff_size_EV_bead+1;
        }
    }

//    /*
//    The last elements in the interval for which there are potential neighbour skips.
//    */
    for (int a=a2+1;a<=A2;a++) {
        int ns=neighbour_skip;
        if (B1<A2) {
            b1 = B1+num_EV;
            if (additional_boundcheck) ns++;
        }
        else {
            b1 = B1;
        }

        int diff = b1-a;
        if (diff <= ns) {
            b1 = a+ns+1;
        }
        b = pmod(b1,num_EV);
        while (b<=B2) {
            dist = singleMove(EV_beads(a),EV_beads(b));
            if (dist < EV_dist) {
                return false;
            }
//            b++;
            b += (dist-EV_dist)/eff_size_EV_bead+1;
        }
    }
    return true;
}
#endif

///////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////
////// CHECK INTERVAL DOUBLE MOVE


#ifndef EXVOL_USE_OLD_CHECK_INTERVAL
bool ExVol::check_interval_doubleMove(int A1, int A2, int B1, int B2) {
/*
    Checks for projected overlap between EV beads in the two intervals limited by
    A1,A2 and B1,B2 respectively. This method assumes that the beads in the B
    interval have not moved, such that only the linear displacement of the beads
    in the A interval have to be considered. The checking is done by the
    doubleMove method.
*/

    #ifdef DEBUG_EXVOL
    if (!(A2<B1 || B2 < A1)) {
        std::cout << "Check Intervals " << A1 << " " << A2 << " " << B1 << " " << B2 << std::endl;
        throw std::invalid_argument("ExVol::check_interval_singleMove(): Invalid Intervals! The invervals are not allowed to cross the periodic boundary!");
    }
    #endif

    double dist;
    int a1,a2;

    if (closed_topology) {
    /*
    ///////////////////////////////////////////////////////////////////////////////////////////////
        PERIODIC BOUNDARY CONDITION
    */

        /*
            Left to right order: A - B
        */
        if (A2<B1) {
            a1 = larger(A1,B2-num_EV+neighbour_skip_boundary_plus_one);
            a2 = smaller(A2,B1-neighbour_skip_plus_one);
            /*
                Check for overlaps within the boundary region of A and region B. Pairs within neighbour skip region are omitted.
            */
            int a = A1;
            while (a<=A2) {
                if ( a1<=a && a<=a2) {
                    a=a2+1;
                    continue;
                }
                int b  = larger (B1,a+neighbour_skip_plus_one);
                int b2 = smaller(B2,a+num_EV-neighbour_skip_boundary_plus_one);
                while (b<=b2) {
                    dist = doubleMove(EV_beads(a),EV_beads(b));
                    if (dist < EV_dist) {
                        return false;
                    }
                    //b++;
                    b += (dist-EV_dist)/eff_size_EV_bead+1;
                }
                a++;
            }
        }
        /*
            Left to right order: B - A
        */
        else {
            a1 = larger(A1,B2+neighbour_skip_plus_one);
            a2 = smaller(A2,B1+num_EV-neighbour_skip_boundary_plus_one);
            /*
                Check for overlaps within the boundary region of A and region B. Pairs within neighbour skip region are omitted.
            */
            int a = A1;
            while (a<=A2) {
                if ( a1<=a && a<=a2) {
                    a=a2+1;
                    continue;
                }
                int b  = larger (B1,a-num_EV+neighbour_skip_boundary_plus_one);
                int b2 = smaller(B2,a-neighbour_skip_plus_one);
                while (b<=b2) {
                    dist = doubleMove(EV_beads(a),EV_beads(b));
                    if (dist < EV_dist) {
                        return false;
                    }
                    //b++;
                    b += (dist-EV_dist)/eff_size_EV_bead+1;
                }
                a++;
            }
        }
    }
    else {
    /*
    ///////////////////////////////////////////////////////////////////////////////////////////////
        NON-PERIODIC BOUNDARY CONDITION
    */
        /*
            Left to right order: A - B
        */
        if (A2<B1) {
            /*
                Define the boundary regions in which pairchecking potentially has to be omitted due to proximity along the chain.
                In this case this region is only at the right boundary of the interval A.
            */
            a1 = A1;
            a2 = smaller(A2,B1-neighbour_skip_plus_one);
            /*
                Check for overlaps within this boundary region of A and region B. Pairs within neighbour skip region are omitted
            */
            for (int a=a2+1;a<=A2;a++) {
                int b = larger(B1,a+neighbour_skip_plus_one);
                while (b<=B2) {
                    dist = doubleMove(EV_beads(a),EV_beads(b));
                    if (dist < EV_dist) {
                        return false;
                    }
                    //b++;
                    b += (dist-EV_dist)/eff_size_EV_bead+1;
                }
            }

        }
        /*
            Left to right order: B - A
        */
        else {
            /*
                Define the boundary regions in which pairchecking potentially has to be omitted due to proximity along the chain.
                In this case this region is only at the left boundary of the interval A.
            */
            a1 = larger(A1,B2+neighbour_skip_plus_one);
            a2 = A2;
            /*
                Check for overlaps within this boundary region of A and region B. Pairs within neighbour skip region are omitted
            */
            for (int a=A1;a<a1;a++) {
                int b  = B1;
                int b2 = smaller(B2,a-neighbour_skip_plus_one);
                while (b<=b2) {
                    dist = doubleMove(EV_beads(a),EV_beads(b));
                    if (dist < EV_dist) {
                        return false;
                    }
                    //b++;
                    b += (dist-EV_dist)/eff_size_EV_bead+1;
                }
            }

        }
    }
    /*
    ///////////////////////////////////////////////////////////////////////////////////////////////
        Finally check all pairs outside the range of potential neighbour skips. I.e. the bulk of the intervals.
    */
    for (int a=a1;a<=a2;a++) {
        int b = B1;
        while (b<=B2) {
            dist = doubleMove(EV_beads(a),EV_beads(b));
            if (dist < EV_dist) {
                return false;
            }
            //b++;
            b += (dist-EV_dist)/eff_size_EV_bead+1;
        }
    }
    return true;
}
#endif

#ifdef EXVOL_USE_OLD_CHECK_INTERVAL
bool ExVol::check_interval_doubleMove(int A1, int A2, int B1, int B2) {
/*
    Checks for projected overlap between EV beads in the two intervals limited by
    A1,A2 and B1,B2 respectively. This method should be used when the beads in both
    intervals have potentially moved such that linear displacements of all beads
    need to be considered. The checking is done by the doubleMove method.
*/

    #ifdef DEBUG_EXVOL
    if (!(A2<B1 || B2 < A1)) {
        std::cout << "Check Intervals " << A1 << " " << A2 << " " << B1 << " " << B2 << std::endl;
        assert(A1<=A2 && B1<=B2);
    }
    if (!(A2<B1 || B2 < A1)) {
        std::cout << "Check Intervals " << A1 << " " << A2 << " " << B1 << " " << B2 << std::endl;
        assert(A2<B1 || B2 < A1);
    }
    #endif

    double dist;
    int b,b1,b2,a1,a2;
    int diff;

    /*
    All pairs outside the range of potential neighbour skips.
    */

    a1 = A1+neighbour_skip;
    if (A1==0 && additional_boundcheck) {
        a1++;
    }
    a2 = A2-neighbour_skip;
    if (A2==num_EV-1 && additional_boundcheck) {
        a2--;
    }

    for (int a=a1;a<=a2;a++) {
        b = B1;
        while (b<=B2) {
            dist = doubleMove(EV_beads(a),EV_beads(b));
            if (dist < EV_dist) {
                return false;
            }
//            b++;
            b += (dist-EV_dist)/eff_size_EV_bead+1;
        }
    }

    /*
    The first elements in the interval for which there are potential neighbour skips.
    */
    for (int a=A1;a<a1;a++) {
        int ns=neighbour_skip;
        if (B2>A1) {
            b2 = B2-num_EV;
            if (additional_boundcheck) ns++;
        }
        else {
            b2 = B2;
        }
        int diff = a-b2;
        if (diff <= ns) {
            b2 = a-ns-1;
        }
        b2 = pmod(b2,num_EV);
        b = B1;
        while (b<=b2) {
            dist = doubleMove(EV_beads(a),EV_beads(b));
            if (dist < EV_dist) {
                return false;
            }
//            b++;
            b += (dist-EV_dist)/eff_size_EV_bead+1;
        }
    }

    /*
    The last elements in the interval for which there are potential neighbour skips.
    */
    for (int a=a2+1;a<=A2;a++) {
        int ns=neighbour_skip;
        if (B1<A2) {
            b1 = B1+num_EV;
            if (additional_boundcheck) ns++;
        }
        else {
            b1 = B1;
        }

        int diff = b1-a;
        if (diff <= ns) {
            b1 = a+ns+1;
        }
        b = pmod(b1,num_EV);;
        while (b<=B2) {
            dist = doubleMove(EV_beads(a),EV_beads(b));
            if (dist < EV_dist) {
                return false;
            }
//            b++;
            b += (dist-EV_dist)/eff_size_EV_bead+1;
        }
    }
    return true;
}
#endif



#ifndef EXVOL_USE_OLD_CHECK_WITHIN_INTERVAL
bool ExVol::check_within_interval(int A1, int A2) {
/*
    Checks whether Excluded Volumes are violated within a single interval
    WARNING: This method has not been fully tested yet.
*/
    int b,b2;
    double dist;
    for (int a=A1;a<A2-neighbour_skip;a++) {
        b  = a+neighbour_skip_plus_one;
        b2 = smaller(a-neighbour_skip_boundary_plus_one+num_EV,A2);
        while (b<=b2) {
            dist = doubleMove(EV_beads(a),EV_beads(b));
            if (dist < EV_dist) {
                return false;
            }
            //b++;
            b += (dist-EV_dist)/eff_size_EV_bead+1;
        }
    }
    return true;
}
#endif

#ifdef EXVOL_USE_OLD_CHECK_WITHIN_INTERVAL
bool ExVol::check_within_interval(int A1, int A2) {
/*
    Checks whether Excluded Volumes are violated within a single interval
    WARNING: This method has not been fully tested yet.
*/
    int b,a1,a2;
    double dist;

    int ns = neighbour_skip;
    if (additional_boundcheck) ns++;

    a1 = A2+ns+1-num_EV;
    bool full_rge_check=true;
    if (a1<=A1) {
        full_rge_check=false;
        a1=A1;
    }

    for (int a=a1;a<=A2;a++) {
        b = a+neighbour_skip+1;
        while (b<=A2) {
            dist = doubleMove(EV_beads(a),EV_beads(b));
            if (dist < EV_dist) {
                return false;
            }
            //b++;
            b += (dist-EV_dist)/eff_size_EV_bead+1;
        }
    }

    if (full_rge_check) {
        for (int a=A1;a<a1;a++) {
            a2 = a-ns-1+num_EV;
            if (A2<a2) {
                a2=A2;
            }
            b = a+neighbour_skip+1;
            while (b<=a2) {
                dist = doubleMove(EV_beads(a),EV_beads(b));
                if (dist < EV_dist) {
                    return false;
                }
                //b++;
                b += (dist-EV_dist)/eff_size_EV_bead+1;
            }
        }
    }
    return true;
}
#endif


double ExVol::doubleMove(int id1, int id2) {

    arma::colvec p1,p1p,p2,p2p,Delta_p,Delta_v,dvec;
    double lambda,dist,dist_primes;

    p1  = bp_pos_backup.col(id1);
    p1p = bp_pos   ->   col(id1);
    p2  = bp_pos_backup.col(id2);
    p2p = bp_pos   ->   col(id2);

    dist_primes = arma::norm(p2p-p1p);
    if (dist_primes<EV_dist) {
        return dist_primes;
    }
    Delta_p = p1-p2;
    Delta_v = p1p-p2p-Delta_p;

    // check if both were translated by the same vector
    dist = arma::dot(Delta_v,Delta_v);
    if (dist < 1e-10) {
        return dist_primes;
    }
    lambda = -arma::dot(Delta_p,Delta_v)/dist;
    // check if the closest approach is outside the moved interval
    if (lambda < 0) {
        return arma::norm(Delta_p);
    }
    if (lambda > 1) {
        return dist_primes;
    }
    dvec = Delta_p+lambda*Delta_v;
    dist = arma::norm(dvec);
    return dist;
}


double ExVol::singleMove(int id1, int id2) {

    arma::colvec p1,p1p,p2,v,w1,w2;
    double n_com, nv;

    p1  = bp_pos_backup.col(id1);
    p1p = bp_pos   ->   col(id1);
    p2  = bp_pos   ->   col(id2);

    v = p1p-p1;
    nv = arma::norm(v);
    w1 = p2-p1;
    if (nv<1e-12) {
        return arma::norm(p2-p1);
    }
    v = v/nv;
    w2 = p1p-p2;

    n_com = arma::dot(v,w1);
    if (n_com < 0) {
        return arma::norm(p2-p1);
    }

    if (arma::dot(v,w2) < 0) {
        return arma::norm(p2-p1p);
    }
    return arma::norm( w1-n_com*v );
}



/////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////
////////////////////// Repulsion Plane //////////////////////////////////////////////////////

void ExVol::set_repulsion_plane(bool active) {
    if (chain->topology_closed()) {
        repulsion_plane_active=false;
    }
    else {
        repulsion_plane_active=active;
        if (active) {
            first_last = REPULSIONPLANE_FIRST_AND_LAST;
        }
    }
}

void ExVol::set_repulsion_plane(bool plane_at_first, bool plane_at_last) {
    if (chain->topology_closed()) {
        repulsion_plane_active=false;
    }
    else {
        if (plane_at_first) {
            repulsion_plane_active=true;
            if (plane_at_last) {
                first_last = REPULSIONPLANE_FIRST_AND_LAST;
            }
            else {
                first_last = REPULSIONPLANE_ONLY_FIRST;
            }
        }
        else {
            if (plane_at_last) {
                repulsion_plane_active=true;
                first_last = REPULSIONPLANE_ONLY_LAST;
            }
            else {
                repulsion_plane_active=false;
            }
        }
    }
}

bool ExVol::RP_check_intervals( const std::vector<arma::ivec>& EV_typeB,
                                const std::vector<arma::ivec>& EV_typeC,
                                const std::vector<arma::ivec>& EV_typeD,
                                const std::vector<arma::ivec>& EV_typeE) {

    /*
        This function does not work because typeA intervals can still slip behind
        the last bead. RP requires all beads to be checked.
    */

    // Check TypeB
    for (unsigned tB=0;tB<EV_typeB.size();tB++) {
        if ( !RP_interval(EV_typeB[tB](EV_FROM),EV_typeB[tB](EV_TO)) ) {
            return false;
        }
    }

    // Check TypeC
    for (unsigned tC=0;tC<EV_typeC.size();tC++) {
        if ( !RP_interval(EV_typeC[tC](EV_FROM),EV_typeC[tC](EV_TO)) ) {
            return false;
        }
    }

    // Check TypeD
    for (unsigned tD=0;tD<EV_typeD.size();tD++) {
        if ( !RP_interval(EV_typeD[tD](EV_FROM),EV_typeD[tD](EV_TO)) ) {
            return false;
        }
    }

    // Check TypeE
    for (unsigned tE=0;tE<EV_typeE.size();tE++) {
        if ( !RP_interval(EV_typeE[tE](EV_FROM),EV_typeB[tE](EV_TO)) ) {
            return false;
        }
    }
    return true;
}


bool ExVol::RP_interval(int A, int B) {

    double z_first, z_last,current_z;
    double dist1,dist2;
    int b;

    z_first = bp_pos->col(0)(2);
    z_last  = bp_pos->col(num_bp-1)(2);

    if (A==0) b=1;
    else      b=A;
    if (B==num_EV-1) B = num_EV-2;

    while (b<=B) {
        current_z = bp_pos->col(EV_beads(b))(2);
        dist1 = current_z - z_first;
        dist2 = z_last    - current_z;

        if (dist1 < 0 || dist2 < 0) {
            return false;
        }
        if (dist2<dist1) {
            dist1 = dist2;
        }
        b += (dist1)/eff_size_EV_bead+1;
    }
    return true;
}

bool ExVol::RP_interval_first(int A, int B) {

    double z_first, current_z;
    double dist;
    int b;

    z_first = bp_pos->col(0)(2);

    if (A==0) b=1;
    else      b=A;

    while (b<=B) {
        current_z = bp_pos->col(EV_beads(b))(2);
        dist = current_z - z_first;

        if (dist < 0) {
            return false;
        }
        b += (dist)/eff_size_EV_bead+1;
    }
    return true;
}

bool ExVol::RP_interval_last(int A, int B) {

    double z_last,current_z;
    double dist;
    int b;

    z_last  = bp_pos->col(num_bp-1)(2);
    b=A;
    if (B==num_EV-1) B = num_EV-2;

    while (b<=B) {
        current_z = bp_pos->col(EV_beads(b))(2);
        dist      = z_last    - current_z;

        if (dist < 0) {
            return false;
        }
        b += (dist)/eff_size_EV_bead+1;
    }
    return true;
}


/////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////
////////////////////// TESTING FUNCTIONS ////////////////////////////////////////////////////


bool ExVol::check_moved_interval_order(const std::vector<arma::ivec>* moved) {
    for (unsigned i=1;i<moved->size();i++) {
        if ((*moved)[i](2) >= 0) {
            if ((*moved)[i](0) <= (*moved)[i-1](1) ) {
                std::cout << "Order violation!" << std::endl;
                std::cout << (*moved)[i-1];
                std::cout << (*moved)[i];
            }
        }
    }
    return true;
}


bool ExVol::check_overlap() {
    arma::colvec p1,p2;
    bool overlap = false;

    double dist;
    int b,b2,a1;
    int diff;

    /*
    All pairs outside the range of potential neighbour skips.
    */

    a1 = neighbour_skip;
    if (additional_boundcheck) {
        a1++;
    }
    for (int a=a1;a<num_EV;a++) {
        b = a+neighbour_skip+1;
        while (b<num_EV) {
            p1 = bp_pos->col(EV_beads(a));
            p2 = bp_pos->col(EV_beads(b));

            dist = arma::norm(p1-p2);
            if (dist < EV_dist) {
                std::cout << "Overlap: " << dist << "(" << EV_dist << ")  -> " << a << " " << b << " - " << EV_beads(a) << " " << EV_beads(b) << std::endl;
                overlap = true;
            }
            b++;
        }
    }

    /*
    The first elements for which there are neighbour skips.
    */
    int ns=neighbour_skip;
    if (additional_boundcheck) ns++;

    for (int a=0;a<a1;a++) {
        b2 = pmod(a-ns-1,num_EV);
        b = a+neighbour_skip+1;
        while (b<=b2) {
            p1 = bp_pos->col(EV_beads(a));
            p2 = bp_pos->col(EV_beads(b));

            dist = arma::norm(p1-p2);
            if (dist < EV_dist) {
                std::cout << "Overlap: " << dist << "(" << EV_dist << ")  -> " << a << " " << b << " - " << EV_beads(a) << " " << EV_beads(b) << std::endl;
                overlap = true;
            }
            b++;
        }
    }
    if (additional_boundcheck) {
        for (unsigned i=0;i<addboundpairs.size();i++) {
            p1 = bp_pos->col(addboundpairs[i](0));
            p2 = bp_pos->col(addboundpairs[i](1));
            dist = arma::norm(p1-p2);
            if (dist < EV_dist) {
                std::cout << "Overlap: " << dist << "(" << EV_dist << ")  -> " << addboundpairs[i](0) << " " << addboundpairs[i](1)  << std::endl;
                std::cout << "special check" << std::endl;
                overlap = true;
            }
        }
    }
    return overlap;
}

