#include "ExcludedVolume.h"

ExVol::ExVol(Chain * ch, double EV_distance, int EV_per_PC) :
chain(ch),
BPS(*ch->get_BPS()),
bp_pos(ch->get_bp_pos()),
triads(ch->get_triads()),
num_bp(ch->get_num_bp()),
num_bps(ch->get_num_bps()),
disc_len(ch->get_disc_len()),
counter(0),
counter_reject(0),
counter_RP_reject(0),
repulsion_plane_active(false),
EV_dist(EV_distance),
num_EV_per_PC(EV_per_PC)
{
    bp_pos_backup = *chain->get_bp_pos();
    triads_backup = *chain->get_triads();

    ///////////////////////////////////////////
    // Init EV_beads

    // set the number of bp per EV bead to the largest odd number smaller than EV_dist/disc_len
    double rest = fmod(EV_dist/disc_len,2);
    if (rest < 1) {
        // even integer
        num_bp_per_EV = int( ( (double)EV_dist/disc_len ) - 1 );
    }
    else {
        // odd integer
        num_bp_per_EV = int( (double)EV_dist/disc_len );
    }

    eff_size_EV_bead = num_bp_per_EV*disc_len;
    num_EV     = int(ceil((double)num_bp/num_bp_per_EV));

    bpid_first_EV   = (num_bp_per_EV-1)/2;
    EV_beads    = arma::zeros<arma::ivec>(num_EV);
    for (int i=0;i<num_EV;i++) {
        int index = bpid_first_EV+i*num_bp_per_EV;
        if (index >= num_bp) index = num_bp-1;
        EV_beads(i) = index;
    }

    cout << endl;
    cout << "#############################################" << endl;
    cout << "####### INITIALIZING EXCLUDED VOLUMES #######" << endl;
    cout << " Excluded Volume Beads: " << endl;
    cout << "   number of EV beads: " << num_EV << endl;
    cout << "   bp per EV bead:     " << num_bp_per_EV << endl;
    cout << "   Effective size:     " << eff_size_EV_bead << endl;
    cout << "   Exclusion distance: " << EV_dist << endl;

    ///////////////////////////////////////////
    // PC beads

    // set num_EV_per_PC to odd number
    if (num_EV_per_PC%2==0) {
        num_EV_per_PC -= 1;
    }

    num_PC          = ceildiv(num_EV,num_EV_per_PC);
    PC_radius       = 0.5*num_EV_per_PC*eff_size_EV_bead;
    num_bp_per_PC   = num_bp_per_EV*num_EV_per_PC;

    EVid_first_PC = (num_EV_per_PC-1 )/2;
    bpid_first_PC = EV_beads(EVid_first_PC);

    for (int i=0;i<num_PC;i++) {
        int from = i*num_EV_per_PC;
        int to   = (i+1)*num_EV_per_PC-1;
        if  (to >= num_EV) {
            to = num_EV-1;
        }

        int middle = from + (to-from)/2;
        int ref_bp = EV_beads(middle);
        arma::colvec initpos = chain->get_bp_pos()->col(ref_bp);

        ProximityCluster * new_PC = new ProximityCluster(from,to,ref_bp,initpos,PC_radius);
        PC_list.push_back(new_PC);
    }

    cout << " Proximity Cluster: " << endl;
    cout << "   Cluster radius:     " << PC_radius << endl;
    cout << "   number of Clusters: " << num_PC << endl;
    cout << "   EV per Cluster:     " << num_EV_per_PC << endl;
    cout << "   bp per Cluster:     " << num_bp_per_PC << endl;
    cout << endl;

    // Set PC_list in Proximity Clusters
    for (int i=0;i<num_PC;i++) {
        PC_list[i]->set_list(PC_list,i);
    }

    ///////////////////////////
    // Build proximity list

    for (int i=0;i<num_PC;i++) {
        PC_list[i]->clear_proximity_list();
    }

    // Calculate proximities
    for (int i=0;i<num_PC;i++) {
        PC_list[i]->calculate_proximity_ascending();
    }

    // build proximity intervals
    for (int i=0;i<num_PC;i++) {
        PC_list[i]->build_proximity_interval();
    }

    steps_since_last_PC_init = 0;
}

ExVol::~ExVol() {

}

double ExVol::get_EV_dist() {
    return EV_dist;
}



inline int ExVol::bp2EV(int bp_id) {
/*
    returns the id of the EV bead that contains the bp with this bp_id
*/
    return bp_id/num_bp_per_EV;
}

inline int ExVol::bp2EV_upper(int bp_id) {
/*
    returns the id of the first EV bead that has a bp index higher or
    equal bp_id
*/
    return (bp_id+bpid_first_EV)/num_bp_per_EV;
}

inline int ExVol::bp2EV_lower(int bp_id) {
/*
    returns the id of the first EV bead that has a bp index smaller or
    equal bp_id
*/
    int check_bead_id = bp_id-bpid_first_EV;
    if (check_bead_id<0) return (bp_id-bpid_first_EV-num_bp_per_EV)/num_bp_per_EV;
    return check_bead_id/num_bp_per_EV;
}

inline int ExVol::bp2PC(int bp_id) {
/*
    returns the id of the Proximity Cluster (PC) that contains the bp with this bp_id
*/
    return bp_id/num_bp_per_PC;
}

inline int ExVol::bp2PC_upper(int bp_id) {
/*
    returns the id of the first Proximity Cluster that has a bp index higher or
    equal bp_id
*/
    return (bp_id+bpid_first_PC)/num_bp_per_PC;
}

inline int ExVol::bp2PC_lower(int bp_id) {
/*
    returns the id of the first Proximity CLuster that has a bp index smaller or
    equal bp_id
*/
    int check_bead_id = bp_id-bpid_first_PC;
    if (check_bead_id<0) return (bp_id-bpid_first_PC-num_bp_per_PC)/num_bp_per_PC;
    return check_bead_id/num_bp_per_PC;
}

inline int ExVol::EV2PC(int EV_id) {
/*
    returns the id of the Proximity Cluster (PC) that contains the EV with this EV_id
*/
    return EV_id/num_EV_per_PC;
}

inline int ExVol::EV2PC_upper(int EV_id) {
/*
    returns the id of the first Proximity Cluster that has a EV index higher or
    equal EV_id
*/
    return (EV_id+EVid_first_PC)/num_EV_per_PC;
}

inline int ExVol::EV2PC_lower(int EV_id) {
/*
    returns the id of the first Proximity CLuster that has a EV index smaller or
    equal EV_id
*/
    int check_bead_id = EV_id-EVid_first_PC;
    if (check_bead_id<0) return (EV_id-EVid_first_PC-num_EV_per_PC)/num_EV_per_PC;
    return check_bead_id/num_EV_per_PC;
}


arma::ivec ExVol::interval_bp2EV(const arma::ivec& bp_interval) {
/*
    Transforms a bp interval into an EV interval. This interval contains all the EV
    beads for which the reference bp (by which the positions of the EV is defined)
    lies inside the given bp interval.
*/
    return {bp2EV_upper(bp_interval(0)),bp2EV_lower(bp_interval(1)),bp_interval(2)};
}

arma::ivec ExVol::interval_bp2PC(const arma::ivec& bp_interval) {
/*
    TODO: Rethink how these intervals have to be transformed. Maybe different
    conversions are necessary depending on the situational need. For the calculation
    of EV violations the PC will have to be included as soon as just one bp falls within
    the range of this PC. However, for the recalculation of the PC proximity based on
    a moved interval a PC should only be recalculated if the reference bp lies within
    the interval.
*/
    return {bp2PC(bp_interval(0)),bp2PC(bp_interval(1)),bp_interval(2)};
}

arma::ivec ExVol::interval_EV2PC(const arma::ivec& EV_interval) {
/*
    TODO: same remark as for interval_bp2PC.
*/
    return {EV2PC(EV_interval(0)),EV2PC(EV_interval(1)),EV_interval(2)};
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


void ExVol::transfer_config(const vector<arma::ivec>* moved, arma::mat* bp_pos_from, arma::mat* bp_pos_to, arma::cube* triads_from, arma::cube* triads_to ) {
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
    *bp_pos = bp_pos_backup;
    *triads = triads_backup;
}

double ExVol::rejection_rate() {
    if (counter>0) {
        return (double)counter_reject/(double)counter;
    }
    else {
        return 0;
    }
}





bool ExVol::check(const vector<arma::ivec>* moved) {
    /*
        TODO: REPULSION PLANES
        - use ProximityClusters
        - integrate plane proximity to ProximityClusters
    */

    counter++;
//    update_PC(moved);
//    check_moved_interval_order(moved);

    vector<arma::ivec> EV_typeA;
    vector<arma::ivec> EV_typeB;
    vector<arma::ivec> EV_typeC;
    vector<arma::ivec> EV_typeD;
    vector<arma::ivec> EV_typeE;
    cal_EV_intervals(moved,EV_typeA,EV_typeB,EV_typeC,EV_typeD,EV_typeE);

//    bool check = check_full();
    bool check = check_intervals(EV_typeA,EV_typeB,EV_typeC,EV_typeD,EV_typeE);

    if (repulsion_plane_active && check) {
        check = RP_check_intervals(EV_typeB,EV_typeC,EV_typeD,EV_typeE);
//        if (!check) {
//            cout << "Repulsion Plane violated!" << endl;
//        }
    }



#ifdef DEBUG_EXVOL
    if (check) {
        if (check_overlap()) {
            cout << "Bead Overlap!" << endl;
            cout << "Moved Intervals:" << endl;
            cout << "Type A" << endl;
            for (int i=0;i<EV_typeA.size();i++) {
                cout << (EV_typeA)[i].t();
            }
            cout << "Type B" << endl;
            for (int i=0;i<EV_typeB.size();i++) {
                cout << (EV_typeB)[i].t();
            }
            cout << "Type C" << endl;
            for (int i=0;i<EV_typeC.size();i++) {
                cout << (EV_typeC)[i].t();
            }
            cout << "Type D" << endl;
            for (int i=0;i<EV_typeD.size();i++) {
                cout << (EV_typeD)[i].t();
            }
            cout << "Type E" << endl;
            for (int i=0;i<EV_typeE.size();i++) {
                cout << (EV_typeE)[i].t();
            }
            std::exit(1);
        }
    }
#endif

    if (check) {
//        cout << "accepted!" << endl;
        /*
            EV respected. Copy config to backup
        */
        transfer_config(moved, bp_pos, &bp_pos_backup, triads, &triads_backup);
    }
    else {
//        cout << "violated!" << endl;
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
    if (counter%EV_FULLCHECK_VIOLATION == 0) {
        if (check_overlap()) {
            /*
                TODO: Add appropriate action
            */
        }
    }

    return check;
}



void ExVol::cal_EV_intervals(   const vector<arma::ivec>* moved,
                                vector<arma::ivec>& EV_typeA,
                                vector<arma::ivec>& EV_typeB,
                                vector<arma::ivec>& EV_typeC,
                                vector<arma::ivec>& EV_typeD,
                                vector<arma::ivec>& EV_typeE) {
/*
    Converts the bp intervals into EV intervals and corrects for potential overlap
    - It is critical that the intervals in 'moved' are ordered from left to right.
*/

    // conversion
    vector<arma::ivec> EV_intervals;
    for (int i=0;i<moved->size();i++) {
        if ((*moved)[i](EV_TYPE) >= 0 && (*moved)[i](EV_FROM) <= (*moved)[i](EV_TO)) {
            EV_intervals.push_back(interval_bp2EV((*moved)[i]));
        }
    }

    // remove overlap
    for (int i=1;i<EV_intervals.size();i++) {
        if (EV_intervals[i](EV_FROM)<=EV_intervals[i-1](EV_TO)) {
            if ( EV_intervals[i](EV_TYPE) >= EV_TYPE_C_FROM ) {
                EV_intervals[i-1](EV_TO) = EV_intervals[i](EV_FROM)-1;
            }
            else {
                EV_intervals[i](EV_FROM) = EV_intervals[i-1](EV_TO)+1;
            }
        }
    }

//    for (int i=0;i<EV_intervals.size();i++) {
//        cout << EV_intervals[i].t() << endl;
//    }

    // sort intervals
    for (int i=0;i<EV_intervals.size();i++) {
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
            cout << "intervals inconsistent!"  << endl;
            for (int j=0;j<EV_intervals.size();j++) {
                cout << EV_intervals[j].t();
            }
            cout << endl;
            for (int i=0;i<moved->size();i++) {
                cout << (*moved)[i].t();
            }
            break;
        }
    }
    if (EV_intervals[EV_intervals.size()-1](1) != num_EV-1) {
        cout << "intervals inconsistent! (last)" << endl;
        cout << "Index of last EV bead: " << num_EV-1 << endl;
        for (int j=0;j<EV_intervals.size();j++) {
            cout << EV_intervals[j].t();
        }
        cout << endl;
        for (int i=0;i<moved->size();i++) {
            cout << (*moved)[i].t();
        }
    }
    #endif
}




bool ExVol::update_PC(const vector<arma::ivec>* moved) {

    arma::ivec   PC_interval;
    arma::colvec new_pos;
    int  last  = -1;
    int  first;

    for (int i=0;i<num_PC;i++) {
        new_pos = chain->get_bp_pos()->col(PC_list[i]->get_ref_bp());
        PC_list[i]->set_pos(new_pos);
    }
    for (int i=0;i<num_PC;i++) {
        PC_list[i]->clear_proximity_list();
    }
    for (int i=0;i<num_PC;i++) {
        PC_list[i]->calculate_proximity_ascending();
    }
    for (int i=0;i<num_PC;i++) {
        PC_list[i]->build_proximity_interval();
    }
    steps_since_last_PC_init=0;
    return false;

    steps_since_last_PC_init++;
    // verify if reinitialization of the positions is necessary
    for (int k=0;k<moved->size();k++) {
        if ((*moved)[k](2) > 0 ) {
            PC_interval = interval_bp2PC((*moved)[k]);
            first = larger(PC_interval(0),last+1);
            last  = PC_interval(1);

            for (int l=first;l<=last;l++) {
                new_pos = chain->get_bp_pos()->col(PC_list[l]->get_ref_bp());
                if ( !PC_list[l]->check_displacement(new_pos) ) {
//                    cout << "Exceeded Max Move!" << endl;

                    for (int i=0;i<num_PC;i++) {
                        new_pos = chain->get_bp_pos()->col(PC_list[i]->get_ref_bp());
                        PC_list[i]->set_pos(new_pos);
                    }
                    for (int i=0;i<num_PC;i++) {
                        PC_list[i]->clear_proximity_list();
                    }
                    for (int i=0;i<num_PC;i++) {
                        PC_list[i]->calculate_proximity_ascending();
                    }
                    for (int i=0;i<num_PC;i++) {
                        PC_list[i]->build_proximity_interval();
                    }
                    steps_since_last_PC_init=0;
                    return false;
                }
            }
        }
    }
    return true;
}


bool ExVol::check_intervals(const vector<arma::ivec>& EV_typeA,
                            const vector<arma::ivec>& EV_typeB,
                            const vector<arma::ivec>& EV_typeC,
                            const vector<arma::ivec>& EV_typeD,
                            const vector<arma::ivec>& EV_typeE) {

    // Check EV_typeC
    for (int tC=0;tC<EV_typeC.size();tC++) {
        // Check with all EV_typeA with singleMove check
        for (int tA=0;tA<EV_typeA.size();tA++) {
            if (!check_singleMove(EV_typeC[tC](EV_FROM),EV_typeC[tC](EV_TO),EV_typeA[tA](EV_FROM),EV_typeA[tA](EV_TO))) {
                return false;
            }
        }
        // Check with all EV_typeB with doubleMove check
        for (int tB=0;tB<EV_typeB.size();tB++) {
            if (!check_doubleMove(EV_typeC[tC](EV_FROM),EV_typeC[tC](EV_TO),EV_typeB[tB](EV_FROM),EV_typeB[tB](EV_TO))) {
                return false;
            }
        }
        // Check with all EV_typeD with doubleMove check
        for (int tD=0;tD<EV_typeD.size();tD++) {
            if (!check_doubleMove(EV_typeC[tC](EV_FROM),EV_typeC[tC](EV_TO),EV_typeD[tD](EV_FROM),EV_typeD[tD](EV_TO))) {
                return false;
            }
        }
        // Check with all EV_typeE with doubleMove check
        for (int tE=0;tE<EV_typeE.size();tE++) {
            if (!check_doubleMove(EV_typeC[tC](EV_FROM),EV_typeC[tC](EV_TO),EV_typeE[tE](EV_FROM),EV_typeE[tE](EV_TO))) {
                return false;
            }
        }
        // Check with all other EV_typeC
        for (int tC2=tC+1;tC2<EV_typeC.size();tC2++) {
            /*
                A single interval can be split into two intervals when it crosses the periodic boundary. For typeC intervals
                these intervals should not be mutually checked
            */
            if ( EV_typeC[tC](EV_TYPE) != EV_typeC[tC2](EV_TYPE) ) {
                if (!check_doubleMove(EV_typeC[tC](EV_FROM),EV_typeC[tC](EV_TO),EV_typeC[tC2](EV_FROM),EV_typeC[tC2](EV_TO))) {
                    return false;
                }
            }
        }
    }

    // Check EV_typeD
    for (int tD=0;tD<EV_typeD.size();tD++) {
        // Check with all EV_typeA with singleMove check
        for (int tA=0;tA<EV_typeA.size();tA++) {
            if (!check_singleMove(EV_typeD[tD](EV_FROM),EV_typeD[tD](EV_TO),EV_typeA[tA](EV_FROM),EV_typeA[tA](EV_TO))) {
//                cout << "D-A violation!" << endl;
//                cout << EV_typeD[tD].t();
//                cout << EV_typeA[tA].t();
                return false;
            }
        }
        // Check with all EV_typeB with doubleMove check
        for (int tB=0;tB<EV_typeB.size();tB++) {
            if (!check_doubleMove(EV_typeD[tD](EV_FROM),EV_typeD[tD](EV_TO),EV_typeB[tB](EV_FROM),EV_typeB[tB](EV_TO))) {
//                cout << "D-B violation!" << endl;
                return false;
            }
        }
        // Check with all EV_typeE with doubleMove check
        for (int tE=0;tE<EV_typeE.size();tE++) {
            if (!check_doubleMove(EV_typeD[tD](EV_FROM),EV_typeD[tD](EV_TO),EV_typeE[tE](EV_FROM),EV_typeE[tE](EV_TO))) {
//                cout << "D-E violation!" << endl;
//                cout << EV_typeD[tD].t();
//                cout << EV_typeE[tE].t();
                return false;
            }
        }

        // Check with all other EV_typeD
        for (int tD2=tD+1;tD2<EV_typeD.size();tD2++) {
            /*
                A single interval can be split into two intervals when it crosses the periodic boundary. For typeD intervals
                these intervals SHOULD be mutally checked regardless!
            */
            if (!check_doubleMove(EV_typeD[tD](EV_FROM),EV_typeD[tD](EV_TO),EV_typeD[tD2](EV_FROM),EV_typeD[tD2](EV_TO))) {
//                cout << "D-D violation!" << endl;
                return false;
            }
        }

        // Check within the interval itself
        if (!check_within_interval(EV_typeD[tD](EV_FROM),EV_typeD[tD](EV_TO))) {
//            cout << "Internal D violation!" << endl;
            return false;
        }
    }

    // Check EV_typeB
    for (int tB=0;tB<EV_typeB.size();tB++) {
        // Check with all EV_typeA with singleMove check
        for (int tA=0;tA<EV_typeA.size();tA++) {
            if (!check_singleMove(EV_typeB[tB](EV_FROM),EV_typeB[tB](EV_TO),EV_typeA[tA](EV_FROM),EV_typeA[tA](EV_TO))) {
                return false;
            }
        }
        // Check with all EV_typeE with doubleMove check
        for (int tE=0;tE<EV_typeE.size();tE++) {
            if (!check_doubleMove(EV_typeB[tB](EV_FROM),EV_typeB[tB](EV_TO),EV_typeE[tE](EV_FROM),EV_typeE[tE](EV_TO))) {
                return false;
            }
        }
    }

    return true;
}


bool ExVol::check_singleMove(int A1, int A2, int B1, int B2) {

//    cout << "Check Intervals " << A1 << " " << A2 << " " << B1 << " " << B2 << endl;
    assert(A2<B1 || B2 < A1);
    assert(A1<=A2 && B1<=B2);
    double dist;
    int b;

    bool FN = false;
    bool BN = false;

    // Forwards Neighbors
    if ( (A2+1)%num_EV == B1 ) {
        FN = true;
    }
    // Backwards Neighbors
    if ( (B2+1)%num_EV == A1 ) {
        BN = true;
    }

    // if interval A consists of a single bead
    if (A1==A2) {
        if (FN) b=B1+1;
        else    b=B1;
        if (BN) B2=B2-1;

        while (b<=B2) {
            dist = singleMove(A1,b);
            if (dist < EV_dist) {
                return false;
            }
//            b++;
//            b += dist/eff_size_EV_bead;
//            b += dist/EV_dist;
            b += (dist-EV_dist)/eff_size_EV_bead+1;
        }
    }
    else {
        // compare all middle ones
        for (int a=A1+1;a<=A2-1;a++) {
            b = B1;
            while (b<=B2) {
                dist = singleMove(a,b);
                if (dist < EV_dist) {
                    return false;
                }
//                b++;
//                b += dist/eff_size_EV_bead;
//                b += dist/EV_dist;
                b += (dist-EV_dist)/eff_size_EV_bead+1;
            }
        }

        if (FN) b=B1+1;
        else    b=B1;
        while (b<=B2) {
            dist = singleMove(A2,b);
            if (dist < EV_dist) {
                return false;
            }
//            b++;
//            b += dist/eff_size_EV_bead;
//            b += dist/EV_dist;
            b += (dist-EV_dist)/eff_size_EV_bead+1;
        }

        if (BN) B2=B2-1;
        b=B1;
        while (b<=B2) {
            dist = singleMove(A1,b);
            if (dist < EV_dist) {
                return false;
            }
//            b++;
//            b += dist/eff_size_EV_bead;
//            b += dist/EV_dist;
            b += (dist-EV_dist)/eff_size_EV_bead+1;
        }
    }
    return true;
}


bool ExVol::check_doubleMove(int A1, int A2, int B1, int B2) {

    if (!(A2<B1 || B2 < A1)) {
        cout << "Check Intervals " << A1 << " " << A2 << " " << B1 << " " << B2 << endl;
        assert(A1<=A2 && B1<=B2);
    }

    if (!(A2<B1 || B2 < A1)) {
        cout << "Check Intervals " << A1 << " " << A2 << " " << B1 << " " << B2 << endl;
        assert(A2<B1 || B2 < A1);
    }

    double dist;
    int b;

    bool FN = false;
    bool BN = false;

    // Forwards Neighbors
    if ( (A2+1)%num_EV == B1 ) {
        FN = true;
    }
    // Backwards Neighbors
    if ( (B2+1)%num_EV == A1 ) {
        BN = true;
    }

    // if interval A consists of a single bead
    if (A1==A2) {
        if (FN) b=B1+1;
        else    b=B1;
        if (BN) B2=B2-1;

        while (b<=B2) {
            dist = doubleMove(A1,b);
            if (dist < EV_dist) {
                return false;
            }
//            b++;
//            b += dist/eff_size_EV_bead;
//            b += dist/EV_dist;
            b += (dist-EV_dist)/eff_size_EV_bead+1;
        }
    }
    else {
        // compare all middle ones
        for (int a=A1+1;a<=A2-1;a++) {
            b = B1;
            while (b<=B2) {
                dist = doubleMove(a,b);
                if (dist < EV_dist) {
                    return false;
                }
//            b++;
//            b += dist/eff_size_EV_bead;
//            b += dist/EV_dist;
            b += (dist-EV_dist)/eff_size_EV_bead+1;
            }
        }

        if (FN) b=B1+1;
        else    b=B1;
        while (b<=B2) {
            dist = doubleMove(A2,b);
            if (dist < EV_dist) {
                return false;
            }
//            b++;
//            b += dist/eff_size_EV_bead;
//            b += dist/EV_dist;
            b += (dist-EV_dist)/eff_size_EV_bead+1;
        }

        if (BN) B2=B2-1;
        b=B1;
        while (b<=B2) {
            dist = doubleMove(A1,b);
            if (dist < EV_dist) {
                return false;
            }
//            b++;
//            b += dist/eff_size_EV_bead;
//            b += dist/EV_dist;
            b += (dist-EV_dist)/eff_size_EV_bead+1;
        }
    }
    return true;
}

bool ExVol::check_within_interval(int A1, int A2) {
/*
    Checks whether Excluded Volumes are violated within a single interval
*/
    int b;
    double dist;

    // check everything but the first
    for (int a=A1+1;a<=A2;a++) {
        b = a+2;
        while (b<=A2) {
            dist = doubleMove(a,b);
            if (dist < EV_dist) {
                return false;
            }
//            b++;
            b += (dist-EV_dist)/eff_size_EV_bead+1;
        }
    }

    // check the first
    if (A1==0 && A2==num_bp-1) {
        A2 = A2-1;
    }
    b = A1+2;
    while (b<=A2) {
        dist = doubleMove(A1,b);
        if (dist < EV_dist) {
            return false;
        }
//        b++;
        b += (dist-EV_dist)/eff_size_EV_bead+1;
    }
    return true;
}


double ExVol::check_interval_boundary(int A, int B) {
/*
    Checks if the boundaries of the interval overlap with any EV_bead
    of the interval
*/
    int b;
    double dist;
    // check the first bead
    b = A+2;
    while (b<=B) {
        dist = doubleMove(A,b);
        if (dist < EV_dist) {
            return false;
        }
//        b++;
        b += (dist-EV_dist)/eff_size_EV_bead+1;
    }
    // check the last bead
    b = A;
    while (b<=B-2) {
        dist = doubleMove(b,B);
        if (dist < EV_dist) {
            return false;
        }
//        b++;
//        b += dist/eff_size_EV_bead;
//        b += dist/EV_dist;
        b += (dist-EV_dist)/eff_size_EV_bead+1;
    }
    return true;
}


double ExVol::doubleMove(int id1, int id2) {

    arma::colvec p1,p1p,p2,p2p,Delta_p,Delta_v,dvec;
    double lambda,dist,dist_primes;

    p1  = bp_pos_backup.col(EV_beads(id1));
    p1p = bp_pos   ->   col(EV_beads(id1));
    p2  = bp_pos_backup.col(EV_beads(id2));
    p2p = bp_pos   ->   col(EV_beads(id2));

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

    p1  = bp_pos_backup.col(EV_beads(id1));
    p1p = bp_pos   ->   col(EV_beads(id1));
    p2  = bp_pos   ->   col(EV_beads(id2));

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
    }
}


bool ExVol::RP_check_intervals( const vector<arma::ivec>& EV_typeB,
                                const vector<arma::ivec>& EV_typeC,
                                const vector<arma::ivec>& EV_typeD,
                                const vector<arma::ivec>& EV_typeE) {

    // Check TypeB
    for (int tB=0;tB<EV_typeB.size();tB++) {
        if ( !RP_interval(EV_typeB[tB](EV_FROM),EV_typeB[tB](EV_TO)) ) {
            return false;
        }
    }

    // Check TypeC
    for (int tC=0;tC<EV_typeC.size();tC++) {
        if ( !RP_interval(EV_typeC[tC](EV_FROM),EV_typeC[tC](EV_TO)) ) {
            return false;
        }
    }

    // Check TypeD
    for (int tD=0;tD<EV_typeD.size();tD++) {
        if ( !RP_interval(EV_typeD[tD](EV_FROM),EV_typeD[tD](EV_TO)) ) {
            return false;
        }
    }

    // Check TypeE
    for (int tE=0;tE<EV_typeE.size();tE++) {
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
//            cout << "dist1 = " << dist1 << endl;
//            cout << "dist2 = " << dist2 << endl;
            return false;
        }
        if (dist2<dist1) {
            dist1 = dist2;
        }
        b += (dist1)/eff_size_EV_bead+1;
    }
    return true;
}






/////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////
////////////////////// TESTING FUNCTIONS ////////////////////////////////////////////////////


bool ExVol::check_moved_interval_order(const vector<arma::ivec>* moved) {
    for (int i=1;i<moved->size();i++) {
        if ((*moved)[i](2) >= 0) {
            if ((*moved)[i](0) <= (*moved)[i-1](1) ) {
                cout << "Order violation!" << endl;
                cout << (*moved)[i-1];
                cout << (*moved)[i];
            }
        }
    }
    return true;
}


bool ExVol::check_full() {
    arma::colvec p1,p2;

    double dist;
    int b;
    for (int a=1;a<num_EV-1;a++) {
        b = a+2;
        while (b<num_EV) {

            dist = doubleMove(a,b);

//            p1 = bp_pos->col(EV_beads(a));
//            p2 = bp_pos->col(EV_beads(b));
//            double diff = arma::norm(p1-p2)-dist;
//            if (diff<-0.5) {
//                cout << "ERROR! -> diff = " << diff << endl;
//                cout << arma::norm(p1-p2) << endl;
//            }


            if (dist < EV_dist) {
                return false;
            }
            b++;
//            b += dist/EV_dist;
//            b += dist/eff_size_EV_bead;
        }
    }

    int a=0;
    b = 2;
    while (b<num_EV-1) {

        dist = doubleMove(a,b);

//        p1 = bp_pos->col(EV_beads(a));
//        p2 = bp_pos->col(EV_beads(b));
//        double diff = arma::norm(p1-p2)-dist;
//        if (diff<-0.5) {
//            cout << "ERROR! -> diff = " << diff << endl;
//            cout << arma::norm(p1-p2) << endl;
//        }


        if (dist < EV_dist) {
            return false;
        }
        b++;
//        b += dist/EV_dist;
//        b += dist/eff_size_EV_bead;
    }
    return true;
}


bool ExVol::check_overlap() {
    arma::colvec p1,p2;
    bool overlap = false;

    double dist;
    int b;
    for (int a=1;a<num_EV-1;a++) {
        b = a+2;
        while (b<num_EV) {

            p1 = bp_pos->col(EV_beads(a));
            p2 = bp_pos->col(EV_beads(b));

            dist = arma::norm(p1-p2);
            if (dist < EV_dist) {
                cout << "Overlap: " << dist << "(" << EV_dist << ")  -> " << a << " " << b << " - " << EV_beads(a) << " " << EV_beads(b) << endl;
                overlap = true;
            }
            b++;
        }
    }

    int a=0;
    b = 2;
    while (b<num_EV-1) {

        p1 = bp_pos->col(EV_beads(a));
        p2 = bp_pos->col(EV_beads(b));

        dist = arma::norm(p1-p2);
        if (dist < EV_dist) {
            cout << "Overlap: " << dist << " -> " << a << " " << b << " - " << EV_beads(a) << " " << EV_beads(b) << endl;
            overlap = true;
        }
        b++;
    }
    return overlap;
}










