#include "ProximityCluster.h"

//////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////// BEAD CLUSTER ////////////////////////////////////////////////////////////////////////

ProximityCluster::ProximityCluster(int EV_first, int EV_last, int reference_bp,const arma::colvec& init_pos, double radius) :
first(EV_first),
last(EV_last),
ref_bp(reference_bp),
num_contained(EV_last-EV_first),
pos(init_pos),
radius(radius)
{
    double prop          = 0.5;
    max_displacement     =   (1+prop)*radius;
    interaction_distance = 2*(1+prop)*radius;
}

void ProximityCluster::set_list(std::vector<ProximityCluster*> bc_list, int own_ID) {
    PC_list = bc_list;
    num_PC  = PC_list.size();
    ID      = own_ID;
}
int  ProximityCluster::get_ID() {
    return ID;
}

int ProximityCluster::get_ref_bp() {
    return ref_bp;
}
int ProximityCluster::get_first() {
    return first;
}
int ProximityCluster::get_last() {
    return last;
}

bool ProximityCluster::check_displacement(arma::colvec& new_pos) {
/*
    Checks of the maximum displacement was exceeded. If so a recalculation
    of the proximity map should be initiated.
*/
    if (arma::norm(pos-new_pos)>max_displacement) {
        return false;
    }
    return true;
}

void ProximityCluster::set_pos(arma::colvec& new_pos) {
/*
    This function does not check the displacement, because positions changes are
    only executed when the maximum displacement was exceeded. Changing this gradually
    would create problems or redundant calculations
*/
    pos = new_pos;
}
arma::colvec* ProximityCluster::get_pos() {
    return &pos;
}

void ProximityCluster::clear_proximity_list() {
    proximity_PCptrs.clear();
    proximity_PCids.clear();
}

void ProximityCluster::calculate_proximity_ascending() {
/*
    Calculates the proximity with respect to all PC to the right in the PC_list. If a PC within
    interaction proximity is found this information is passed to the involved PC.
    Looping through PC_list in ascending order allows for all proximities to be linked, while
    only calculating every pair distance once. With this method proximity_PCptrs will
    contain all PC of potential overlap in ascending order (provided that this method is
    called for ALL PC in ascending order.
*/
    proximity_PCptrs.push_back(PC_list[ID]);
    for (int i=ID+1;i<num_PC;i++) {
        if (arma::norm(*PC_list[i]->get_pos()-pos) < interaction_distance) {
            proximity_PCptrs.push_back(PC_list[i]);
            PC_list[i]->add_to_proximity(ID);
        }
    }
}

void ProximityCluster::add_to_proximity(int id) {
/*
    Allows calculate_proximity_ascending to pass the proximity information.
*/
    proximity_PCptrs.push_back(PC_list[id]);
}

void ProximityCluster::calculate_proximity_full() {
/*
    Find all PC in proximity. Only to be used if only this ProximityCluster is to be updated.
*/
    for (int i=0;i<num_PC;i++) {
        if (arma::norm(*PC_list[i]->get_pos()-pos) < interaction_distance) {
            proximity_PCptrs.push_back(PC_list[i]);
        }
    }
}

void ProximityCluster::build_proximity_interval() {
/*
    Casts the information of ProximityClusters in proximity into
    merged intervals, which allows for more efficient looping, due to
    omitting beads that are impossible to be sufficiently close based
    on the distance of the previous bead. If all the PC in proximity
    are sequentially adjacent, this function will construct a single
    interval.
*/
    int from, to, curr_ID, next_ID;
    proximity_interval.clear();
    num_proximity_intervals = 1;

    from    = proximity_PCptrs[0]->get_first();
    curr_ID = proximity_PCptrs[0]->get_ID();

    for (unsigned i=1;i<proximity_PCptrs.size();i++) {
        next_ID = proximity_PCptrs[i]->get_ID();
        if (next_ID==ID){
            proximity_interval_self = num_proximity_intervals-1;
        }
        if (next_ID!=curr_ID+1) {
            proximity_interval.push_back({from,proximity_PCptrs[i-1]->get_last()});
            from = proximity_PCptrs[i]->get_first();
            num_proximity_intervals++;
        }
        curr_ID=next_ID;
    }
    proximity_interval.push_back({from,proximity_PCptrs[proximity_PCptrs.size()-1]->get_last()});
//    assert(num_proximity_intervals==proximity_interval.size());
}

int  ProximityCluster::get_proximity_interval(std::vector<arma::ivec>** proximity_interval_ptr) {
/*
    Retrieves the vector of proximity intervals. The vector has to be allocated outside this method and a pointer
    to the pointer to this vector has to be passed as an argument.
*/
    *proximity_interval_ptr = &proximity_interval;
    return proximity_interval_self;
}

void ProximityCluster::get_proximity_tail(  std::vector<arma::ivec>** proximity_interval_ptr,
                                            arma::colvec& rel_interval,
                                            int left_limit,
                                            int right_limit)
                                            {
/*
    Provides the proximity intervals relevant for tail segments.
    The tail is assumed to be translated by a single common displacement, which means that mutual overlap
    is impossible. Overlap between these segments and type I and type II intervals is assumed to be treated
    somewhere else. This method identifies the intervals of potential overlap with this PC, which are not
    contained in any of the moved intervals.

    For now it is assumed that non-moved intervals are to the left of the first moved interval. left_limit
    is the index of the first element in the first moved interval.

    The relevant interval is the interval of the Proximity cluster that is not located within a moved interval,
    i.e. to the right of right_limit.

*/

    // find the proximity intervals located to the left of left_limit
    std::vector<arma::ivec> tail_interval;
    int A0,A1;
    for (int pi=0;pi<num_proximity_intervals;pi++) {
        A0 = proximity_interval[pi](0);
        A1 = proximity_interval[pi](1);
        if (A0 <  left_limit) {
            tail_interval.push_back({A0,smaller(A1,left_limit-1)});
        }
        if (A1 >= left_limit) {
            break;
        }
    }
    *proximity_interval_ptr = &tail_interval;

    // find the relevant interval
    if (first <= right_limit) {
        rel_interval = {right_limit+1,last};
    }
    else {
        rel_interval = {first,last};
    }
}








