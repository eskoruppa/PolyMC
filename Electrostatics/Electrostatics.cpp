#include "Electrostatics.h"

ElStat::ElStat(Chain * ch, ESPotential * espot, double rho_max, double rho_min, double neighbor_skip_dist) :
chain(ch),
BPS(*ch->get_BPS()),
bp_pos(ch->get_bp_pos()),
triads(ch->get_triads()),
espot(espot),
num_bp(ch->get_num_bp()),
num_bps(ch->get_num_bps()),
disc_len(ch->get_disc_len()),
kT(chain->get_kT()),
rho_max(rho_max),
rho_min(rho_min),
neighbor_skip_dist(neighbor_skip_dist),
pending_check(false)
{
    if (chain->topology_closed() || chain->topology_pseudo_closed()) {
        closed_topology = true;
    }
    else {
        closed_topology = false;
    }

    ///////////////////////////////////////////
    // Init matrices
    if (closed_topology) {
        num_midpoint = num_bp;
    }
    else {
        num_midpoint = num_bp-1;
    }
    half_num_midpoint = num_midpoint/2;

    midpoint_pos        = arma::zeros(3,num_midpoint);
    midpoint_pos_backup = arma::zeros(3,num_midpoint);

    pair_interactions           = arma::zeros(num_midpoint,num_midpoint);
    pair_interactions_backup    = arma::zeros(num_midpoint,num_midpoint);

    if (std::abs(neighbor_skip_dist-disc_len) < 1e-10) {
        neighbor_skip_dist = disc_len;
    }
    neighbor_skip_num = int(std::floor(neighbor_skip_dist/disc_len));

    cal_midpoint_pos();
    eval_full();
    set_current_as_backup();

    counter = 0;

    #ifdef DEBUG_ELSTAT
    check_pair_interactions();
    #endif

    std::cout << std::endl;
    std::cout << "######################################" << std::endl;
    std::cout << "#### INITIALIZING ELECTROSTATICS #####" << std::endl;
    std::cout << "######################################" << std::endl;
    std::cout << "   potenital type:     " << espot->get_type()  << std::endl;
    std::cout << "   midpoint beads:     " << num_midpoint       << std::endl;
    std::cout << "   neighborskip segs   " << neighbor_skip_num  << std::endl;
    std::cout << "   neighborskip dist   " << neighbor_skip_dist << std::endl;
    std::cout << "   closed topology:    " << closed_topology    << std::endl;
    std::cout << "######################################" << std::endl;
    std::cout << std::endl;
}

ElStat::~ElStat() {

}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////

double ElStat::get_current_energy(bool kt_units,bool recal) {
    if (recal) {
        eval_full();
    }
    double E = arma::accu(pair_interactions)*0.5;
    if (kt_units) {
        E /= kT;
    }
    return E;
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void    ElStat::cal_midpoint_pos(const std::vector<arma::ivec> & midpoint_intervals) {

    unsigned from,to,last;
    /*
        all except last one can be done collectively
    */
    for (unsigned i=0;i<midpoint_intervals.size()-1;i++) {
        from = midpoint_intervals[i](EV_FROM);
        to   = midpoint_intervals[i](EV_TO);
        if (within_EV_moved_type(midpoint_intervals[i](EV_TYPE)) && from <= to) {
            for (unsigned j=from;j<=to;j++) {
                midpoint_pos.col(j) = 0.5*(bp_pos->col(j)+bp_pos->col(j+1));
            }
        }
    }
    /*
        IF CLOSED TOPOLOGY
    */
    last = midpoint_intervals.size()-1;
    if (closed_topology) {
        from = midpoint_intervals[last](EV_FROM);
        to   = midpoint_intervals[last](EV_TO);
        if (within_EV_moved_type(midpoint_intervals[last](EV_TYPE)) && from <= to) {
            to--;
            for (unsigned j=from;j<=to;j++) {
                midpoint_pos.col(j) = 0.5*(bp_pos->col(j)+bp_pos->col(j+1));
            }
            /*
                handle last
            */
            #ifdef DEBUG_ELSTAT
            if (to+1 != num_bp-1) {
                throw std::invalid_argument("ElStat::cal_midpoint_pos(): to+1 != num_bp-1 for closed topology!");
            }
            #endif
            midpoint_pos.col(num_midpoint-1) = 0.5*(bp_pos->col(num_midpoint-1)+bp_pos->col(0));
        }

    }
    /*
        IF OPEN TOPOLOGY
    */
    else {
        /*
            if not closed this interval can simply be calculated
        */
        from = midpoint_intervals[last](EV_FROM);
        to   = midpoint_intervals[last](EV_TO);
        if (within_EV_moved_type(midpoint_intervals[last](EV_TYPE)) && from <= to) {
            for (unsigned j=from;j<=to;j++) {
                midpoint_pos.col(j) = 0.5*(bp_pos->col(j)+bp_pos->col(j+1));
            }
        }
    }
}

void    ElStat::cal_midpoint_pos() {
    for (int i=0;i<num_bp-1;i++) {
        midpoint_pos.col(i) = 0.5*(bp_pos->col(i)+bp_pos->col(i+1));
    }
    if (closed_topology) {
        midpoint_pos.col(num_bp-1) = 0.5*(bp_pos->col(num_bp-1)+bp_pos->col(0));
    }
}


void    ElStat::cal_midpoint_intervals( const std::vector<arma::ivec>* moved,
                                        std::vector<arma::ivec>& midpoint_intervals) {

    /*
        TODO: KEEP ONLY 3 types of midpoint intervals:
            - not moved
            - moved without internal checks
            - moved with internal check

        combine neighboring intervals of same type. Like this the element skip can be further leveraged!
    */


    /*
        moved may contain EV_TYPE -1, or have smaller upper bound than lower bound indicating
        a zero size interval
    */
    int eff_first,eff_last,eff_next;
    bool lower_next;
    int from,to;

    /*
        FIND EFFECTIVE FIRST AND LAST MOVED INTERVALS
    */
    eff_first = -1;
    eff_last  = -1;
    for (unsigned i=0;i<moved->size();i++) {
        if (within_EV_active_type((*moved)[i](EV_TYPE)) && (*moved)[i](EV_FROM) <= (*moved)[i](EV_TO)) {
            eff_first = i;
            break;
        }
    }
    for (unsigned i=moved->size()-1;i>=0;i--) {
        if (within_EV_active_type((*moved)[i](EV_TYPE)) && (*moved)[i](EV_FROM) <= (*moved)[i](EV_TO)) {
            eff_last = i;
            break;
        }
    }
    if (eff_first < 0) {
        return;
    }

    /*
        ASSIGN INTERVALS EXCEPT LAST
    */
    lower_next = false;
    /*
        last is considered separately
    */
    for (unsigned i=eff_first;i<=eff_last-1;i++) {
        if (within_EV_active_type((*moved)[i](EV_TYPE)) && (*moved)[i](EV_FROM) <= (*moved)[i](EV_TO)) {

            from = (*moved)[i](EV_FROM);
            to   = (*moved)[i](EV_TO);
            if (lower_next) {
                from--;
            }
            eff_next = i+1;
            while (!within_EV_active_type((*moved)[eff_next](EV_TYPE))) {
                eff_next++;
                #ifdef DEBUG_ELSTAT
                if (eff_next > eff_last) {
                    throw std::invalid_argument("ElStat::cal_midpoint_intervals(): eff_next exceeds eff_last");
                }
                #endif
            }


            if ((*moved)[i](EV_TYPE) < (*moved)[eff_next](EV_TYPE)) {
                lower_next = true;
                to--;
            }
            else{
                lower_next = false;
            }
            if (from<=to) {
                arma::ivec new_interval = {from,to,(*moved)[i](EV_TYPE)};
                midpoint_intervals.push_back(new_interval);
            }
        }
    }

    /*
        ASSIGN LAST INTERVALS
    */
    from = (*moved)[eff_last](EV_FROM);
    to   = (*moved)[eff_last](EV_TO);
    if (lower_next) {
        from--;
    }

    #ifdef DEBUG_ELSTAT
    if (to!=num_bp-1) {
        throw std::invalid_argument("ElStat::cal_midpoint_intervals(): Invalid Intervals! Last inteval does not extend to the last bead");
    }
    #endif

    if (closed_topology) {
        if ((*moved)[eff_last](EV_TYPE) < (*moved)[eff_first](EV_TYPE)) {
            if (from<to) {
                arma::ivec new_interval = {from,to-1,(*moved)[eff_last](EV_TYPE)};
                midpoint_intervals.push_back(new_interval);
            }
            // add additional interval of size 1 and same type as first interval
            arma::ivec new_interval = {to,to,(*moved)[eff_first](EV_TYPE)};
            midpoint_intervals.push_back(new_interval);
        }
        else {
            /*
                if first doesn't overrule last, simply add last over given range
            */
            if (from<=to) {
                arma::ivec new_interval = {from,to,(*moved)[eff_last](EV_TYPE)};
                midpoint_intervals.push_back(new_interval);
            }
        }
    }
    else {
        /*
            If topology is not closed, there are num_bp-1 midpoint beads!
        */
        if (from<to) {
            arma::ivec new_interval = {from,to-1,(*moved)[eff_last](EV_TYPE)};
            midpoint_intervals.push_back(new_interval);
        }
    }

    #ifdef DEBUG_ELSTAT
    if (midpoint_intervals[0](EV_FROM) != 0) {
        throw std::invalid_argument("ElStat::cal_midpoint_intervals(): from of first interval is not zero");
    }
    for (unsigned i=0;i<midpoint_intervals.size()-1;i++) {
        if (midpoint_intervals[i+1](EV_FROM) - midpoint_intervals[i](EV_TO) != 1) {
            throw std::invalid_argument("ElStat::cal_midpoint_intervals(): difference different from one");
        }
    }
    if (midpoint_intervals[midpoint_intervals.size()-1](EV_TO) != num_midpoint-1) {
        throw std::invalid_argument("ElStat::cal_midpoint_intervals(): to of last not equal to num_midpoint-1");
    }
    #endif
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////


double  ElStat::cal_beta_dE(const std::vector<arma::ivec>* moved){

    if (pending_check) {
        std::cout << "Error: ElStat::cal_beta_dE -> last check still pending!" << std::endl;
    }
    pending_check = true;
    counter++;

    std::vector<arma::ivec> midpoint_intervals;
    cal_midpoint_intervals( moved,midpoint_intervals);
    cal_midpoint_pos(midpoint_intervals);

    #ifdef DEBUG_ELSTAT
    check_midpoint_pos();
    #endif

//    std::cout << " #################### " << std::endl;
//    std::cout << " #################### " << std::endl;
//    std::cout << " INTERVALS " << std::endl;
//    for (unsigned i=0;i<midpoint_intervals.size();i++) {
//        std::cout << midpoint_intervals[i](EV_TYPE) << " " << midpoint_intervals[i](EV_FROM) << " " << midpoint_intervals[i](EV_TO) << std::endl;
//    }
//    std::cout << " #################### " << std::endl;


    double dE = 0;
    for (unsigned i=0;i<midpoint_intervals.size();i++) {
        if (within_EV_typeA(midpoint_intervals[i](EV_TYPE))) {
            // only check moved intervals
            for (unsigned j=i+1;j<midpoint_intervals.size();j++) {
                if (!within_EV_typeA(midpoint_intervals[j](EV_TYPE))) {
                    dE += eval_dE_inter_interval(   midpoint_intervals[i](EV_FROM),
                                                    midpoint_intervals[i](EV_TO),
                                                    midpoint_intervals[j](EV_FROM),
                                                    midpoint_intervals[j](EV_TO));
                }
            }
        }
        else {
            // check all
            if (i < midpoint_intervals.size()-1) {
                dE += eval_dE_inter_interval(   midpoint_intervals[i](EV_FROM),
                                                midpoint_intervals[i](EV_TO),
                                                midpoint_intervals[i](EV_TO)+1,
                                                num_midpoint-1);
            }

            if (within_EV_typeD(midpoint_intervals[i](EV_TYPE))) {
                /*
                    check all internally
                */
                dE += eval_dE_intra_interval(   midpoint_intervals[i](EV_FROM),
                                                midpoint_intervals[i](EV_TO));
            }
            else {
                /*
                    check only boundary beads internally
                */
                int a,b;
                a = midpoint_intervals[i](EV_FROM);
                b = midpoint_intervals[i](EV_TO);
                dE += eval_dE_bead_on_interval(a,a+1,b);
                dE += eval_dE_bead_on_interval(b,a+1,b-1);
            }

            #ifdef DEBUG_ELSTAT
            if (within_EV_typeE(midpoint_intervals[i](EV_TYPE))) {
                throw std::invalid_argument("ElStat::calculate_electrostatics(): electrostatics not yet compatible with typeE moves (index translations)");
            }
            #endif

        }
    }

    #ifdef DEBUG_ELSTAT
    check_pair_interactions();
    #endif
    return dE/kT;
}

void  ElStat::eval_full(){
    eval_dE_intra_interval(0, num_midpoint-1);
    set_current_as_backup();
}


double ElStat::eval_dE_inter_interval(int a1, int a2, int b1, int b2) {
    /*
        assumes a1 <= a2 < b1 <= b2
    */
    double dE = 0;
    for (int a=a1;a<=a2;a++) {
        dE += eval_dE_bead_on_interval(a, b1, b2);
    }
    return dE;
}

double ElStat::eval_dE_intra_interval(int a1, int a2) {
    double dE = 0;
    for (int a=a1;a<a2;a++) {
        dE += eval_dE_bead_on_interval(a, a+1, a2);
    }
    return dE;
}

double ElStat::eval_dE_bead_on_interval(int a, int b1, int b2) {
    /*
        assumes a < b1 <= b2
    */

    #ifdef DEBUG_ELSTAT
    if ((a >= b1 && a <= b2)){
        std::cout << a << " " << b1 << " " << b2 << std::endl;
        throw std::invalid_argument("ElStat::eval_bead_on_inteval(): error in interval ordering");
    }
    #endif

    if (b2 < b1) {
        return 0;
    }

    int nb1,nb2;
    #ifdef USE_STEPSKIP
    int bn;
    #endif
    nb1 = impose_skip(a,b1);
    nb2 = impose_skip(a,b2);

    if (nb1<b1) {
        return 0;
    }
    if (nb2>b2) {
        return 0;
    }
    if (nb2<0) {
        return 0;
    }
    if (nb1 >= num_midpoint) {
        return 0;
    }

    double dE = 0;
    double dist,energy;
    int b = nb1;
    while (b<=nb2) {
        dist = arma::norm(midpoint_pos.col(a)-midpoint_pos.col(b));
        if (dist<=rho_max) {
            energy = eval_pair(a, b, dist);

            #ifdef DEBUG_ELSTAT
            double rev = eval_pair(b,a, dist);
            if (std::abs(energy-rev) > 1e-10) {
                if (espot->tabulating()) {
                    if (std::abs(energy-rev) < 1e-1) {
                        continue;
                    }
                }
                std::cout << "rev energy inconsistent" << std::endl;
                std::cout << energy << std::endl;
                std::cout << rev << std::endl;
                std::exit(0);
            }
            #endif
            pair_interactions(a,b) = energy;
            pair_interactions(b,a) = energy;
            dE += energy - pair_interactions_backup(a,b);
            b += 1;
        }
        else {
            #ifdef USE_STEPSKIP
            // USE_STEPSKIP
            bn = b+std::floor((dist-rho_max)/disc_len)+1;
            if (bn>nb2) {
                bn = nb2+1;
            }
            for (int bs=b;bs<bn;bs++) {
                #ifdef DEBUG_ELSTAT
                    if (std::abs(eval_pair(a, bs, dist))>1e-14) {
                        std::cout << "interaction not zero" << std::endl;
                        std::cout << " " << eval_pair(a, bs, dist) << std::endl;
                    }
                #endif
                pair_interactions(a,bs) = 0;
                pair_interactions(bs,a) = 0;
                dE -= pair_interactions_backup(a,bs);
            }
            b = bn;
            #else
            pair_interactions(a,b) = 0;
            pair_interactions(b,a) = 0;
            dE -= pair_interactions_backup(a,b);
            b  += 1;
            #endif
        }
    }
    return dE;
}

double ElStat::eval_pair(int a, int b) {
    #ifdef DEBUG_ELSTAT
    if (in_neighborskip_distance(a,b)) {
        std::cout << "in neighbor skip distance" << std::endl;
        std::cout << a << " - " << b << std::endl;
        std::exit(0);
    }
    #endif
    return espot->segmentpair_energy(midpoint_pos.col(a)-midpoint_pos.col(b),triads->slice(a).col(2),triads->slice(b).col(2));
}

double ElStat::eval_pair(int a, int b, double dist) {
    #ifdef DEBUG_ELSTAT
    if (in_neighborskip_distance(a,b)) {
        std::cout << "in neighbor skip distance" << std::endl;
        std::cout << a << " - " << b << std::endl;
        std::exit(0);
    }
    #endif
    return espot->segmentpair_energy(midpoint_pos.col(a)-midpoint_pos.col(b),triads->slice(a).col(2),triads->slice(b).col(2));
}




///////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////


int ElStat::impose_skip(int a, int b) {
    int diff;
    if (closed_topology) {
        if (a<=b) {
            diff = b-a;
            if (diff<=neighbor_skip_num) {
                return a+neighbor_skip_num+1;
            }
            else {
                diff = num_midpoint-diff;
                if (diff<=neighbor_skip_num) {
                    return a+num_midpoint-neighbor_skip_num-1;
                }
            }
        }
        else {
            diff = a-b;
            if (diff<=neighbor_skip_num) {
                return a-neighbor_skip_num-1;
            }
            else {
                diff = num_midpoint-diff;
                if (diff<=neighbor_skip_num) {
                    return a-num_midpoint+neighbor_skip_num+1;
                }
            }
        }
    }
    else {
        if (a<=b) {
            if (b-a<=neighbor_skip_num) {
                return a+neighbor_skip_num+1;
            }
        }
        else {
            if (a-b<=neighbor_skip_num) {
                return a-neighbor_skip_num-1;
            }
        }
    }
    return b;
}

bool ElStat::in_neighborskip_distance(int a, int b) {
    int diff;
    if (closed_topology) {
        if (a<=b) {
            diff = b-a;
            if (diff<=neighbor_skip_num) {
                return true;
            }
            else {
                diff = num_midpoint-diff;
                if (diff<=neighbor_skip_num) {
                    return true;
                }
            }
        }
        else {
            diff = a-b;
            if (diff<=neighbor_skip_num) {
                return true;
            }
            else {
                diff = num_midpoint-diff;
                if (diff<=neighbor_skip_num) {
                    return true;
                }
            }
        }
    }
    else {
        if (a<=b) {
            if (b-a<=neighbor_skip_num) {
                return true;
            }
        }
        else {
            if (a-b<=neighbor_skip_num) {
                return true;
            }
        }
    }
    return false;
}


///////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////

/*
    Backup Handling
*/

void ElStat::revert_to_backup() {
    midpoint_pos        = midpoint_pos_backup;
    pair_interactions   = pair_interactions_backup;

    #ifdef DEBUG_ELSTAT
    check_midpoint_pos();
    check_pair_interactions();
    #endif
    pending_check = false;
    if (counter%RECAL_FULL_EVERY==0) {
        recal_full_energy();
    }
}

void ElStat::set_current_as_backup(bool recal) {
    if (recal) {
        cal_midpoint_pos();
        eval_full();
    }
    midpoint_pos_backup        = midpoint_pos;
    pair_interactions_backup   = pair_interactions;

    #ifdef DEBUG_ELSTAT
    check_midpoint_pos();
    check_pair_interactions();
    #endif
    pending_check = false;
    if (counter%RECAL_FULL_EVERY==0) {
        recal_full_energy();
    }
}

///////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////

void ElStat::recal_full_energy() {
    #ifdef DEBUG_ELSTAT
    if (pending_check) {
        std::cout << "Cannot run recal_full_energy when check is pending"  << std::endl;
        std::exit(0);
    }
    check_midpoint_pos();
    #endif
    double e,E;
    E = 0;
    for (unsigned a=0;a<num_midpoint-1;a++) {
        for (unsigned b=a+1+neighbor_skip_num;b<impose_skip(a,num_midpoint);b++) {
            e = eval_pair(a, b);
            if (!equal_double(e,pair_interactions(a,b),MAX_ENERGY_DIFF)) {
                double dist = arma::norm(midpoint_pos_backup.col(a)-midpoint_pos_backup.col(b));
                if (equal_double(dist,rho_max),1e-10) {
                    continue;
                }
                std::cout << "Energy inconsistent"  << std::endl;
                std::cout << a << " - " << b << std::endl;
                std::cout << "current:    " << e << std::endl;
                std::cout << "stored:     " << pair_interactions(a,b) << std::endl;
                std::cout << "inverse:    " << eval_pair(b,a) << std::endl;
                std::cout << " diff:      " << pair_interactions(a,b) - e << std::endl;
                std::cout << " dist:      " << dist << std::endl;
                if (!espot->tabulating()) {
                    throw std::invalid_argument("ElStat::recal_full_energy(): Energy inconsistent");
                }
            }
            pair_interactions(a,b) = e;
            pair_interactions(b,a) = e;
            pair_interactions_backup(a,b) = e;
            pair_interactions_backup(b,a) = e;
        }
    }
}



///////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////

/*
    CHECK FUNCTIONS
*/

bool ElStat::check_midpoint_pos() {
    arma::mat  midpoint_pos_check = arma::zeros(3,num_midpoint);
    for (int i=0;i<num_bp-1;i++) {
        midpoint_pos_check.col(i) = 0.5*(bp_pos->col(i)+bp_pos->col(i+1));
    }
    if (closed_topology) {
        midpoint_pos_check.col(num_bp-1) = 0.5*(bp_pos->col(num_bp-1)+bp_pos->col(0));
    }

    for (int i=0;i<num_midpoint;i++) {
        double diff = arma::norm(midpoint_pos_check.col(i)-midpoint_pos.col(i));
        if (diff > 1e-10) {
            std::cout << "midpoint_pos not correctly calculated!" << std::endl;
            std::exit(0);
            return false;
        }
    }
    return true;
}


double ElStat::check_eval_pair(int a, int b) {
    if (in_neighborskip_distance(a,b)) {
        return 0;
    }
    return espot->segmentpair_energy(midpoint_pos.col(a)-midpoint_pos.col(b),triads->slice(a).col(2),triads->slice(b).col(2));
}

bool ElStat::check_pair_interactions() {
    arma::mat pi_check = arma::zeros(num_bp,num_midpoint);
    double E;

    for (unsigned a=0;a<num_midpoint-1;a++) {
        for (unsigned b=a+1;b<num_midpoint;b++) {
            double E1,E2;
            E1 = pair_interactions(a,b);
            E2 = pair_interactions(b,a);
            if (std::abs(E1-E2)>1e-10) {
                std::cout << "energy assymetric" << std::endl;
                std::cout << a << " - " << b << std::endl;
                std::cout << E1 << std::endl;
                std::cout << E2 << std::endl;
                std::exit(0);
            }
        }
    }


    double e;
    E = 0;
    for (unsigned a=0;a<num_midpoint-1;a++) {
        for (unsigned b=a+1;b<num_midpoint;b++) {
            e = check_eval_pair(a, b);
            if (!equal_double(e,pair_interactions(a,b),1e-10)) {
                std::cout << "Energy inconsistent"  << std::endl;
                std::cout << "check_pair_interactions"  << std::endl;
                std::cout << a << " - " << b << std::endl;
                std::cout << e << std::endl;
                std::cout << pair_interactions(a,b) << std::endl;
                std::cout << pair_interactions(a,b)-e << std::endl;
                std::exit(0);
            }
        }
    }

    ///////
    return true;
}



void ElStat::set_backup_pos(arma::mat*  bp_pos_backup,arma::cube* triads_backup) {
    this->bp_pos_backup = bp_pos_backup;
    this->triads_backup = triads_backup;
}
