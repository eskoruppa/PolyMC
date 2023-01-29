#ifndef __EXCLUDED_VOLUME__
#define __EXCLUDED_VOLUME__

#include "../Chain.h"
#include "../SO3Methods.h"
#include "../ExtraFuncs.h"

#include "../Geometry/Geometry.h"
#include "ClosureCrossings.h"


#define _USE_MATH_DEFINES

#include <armadillo> // armadillo
#include <iostream>  // input output
#include <string>    // string
#include <algorithm>
#include <stdexcept> // exceptions
#include <fstream>   // read file
#include <cmath>     // math lib
#include <vector>    // vector
#include <cassert>   // assert

#define EV_TYPE_A 0
#define EV_TYPE_B 1
#define EV_TYPE_C_FROM 1000
#define EV_TYPE_C_TO   1999
#define EV_TYPE_D_FROM 2000
#define EV_TYPE_D_TO   2999
#define EV_TYPE_E 2

#define EV_FROM     0
#define EV_TO       1
#define EV_TYPE     2

#define EV_FULLCHECK_VIOLATION 1000000

#define REPULSIONPLANE_FIRST_AND_LAST 0
#define REPULSIONPLANE_ONLY_FIRST     1
#define REPULSIONPLANE_ONLY_LAST      2

/*
    DEBUG FLAG
*/
//#define DEBUG_EXVOL
//#define EXVOL_USE_OLD_CHECK_INTERVAL
//#define EXVOL_USE_OLD_CHECK_WITHIN_INTERVAL


/*
Different types of moved intervals:


TypeA:
    indicator range:        0
    type of displacement:   Beads that did not move and thus don't require crosschecking with other
                            non-moved beads

TypeB:
    indicator range:        1
    type of displacement:   Beads that were collectively trianslated. Require only single move check with
                            TypeA.

TypeC:
    indicator range:        1000-1999
    type of displacement:   Collectively moved beads that don't require internal checking but require double
                            move checks with other TypeC and TypeD intervals

TypeD:
    indicator range:        2000-2999
    type of displacement:   Intervals in which beads can move individually, requireing crosschecks within the interval
                            itself

TypeE:
    indicator range:        2
    type of displacement:   Interval for beads that are translated onto the former position of another bead. This requires
                            double move checks to prevent false collision detection.

REMARKS:
    The moved intervals passed by the MCStep method have to be in ascending order.
*/

inline bool within_EV_moved_type(int id) {
    return (id > EV_TYPE_A);
}
inline bool within_EV_active_type(int id) {
    return (id >= 0);
}
inline bool within_EV_typeA(int id) {
    return (id == EV_TYPE_A);
}
inline bool within_EV_typeB(int id) {
    return (id == EV_TYPE_B);
}
inline bool within_EV_typeC(int id) {
    return (EV_TYPE_C_FROM <= id && id <= EV_TYPE_C_TO);
}
inline bool within_EV_typeD(int id) {
    return (EV_TYPE_D_FROM <= id && id <= EV_TYPE_D_TO);
}
inline bool within_EV_typeE(int id) {
    return (id == EV_TYPE_E);
}



class ExVol;

class ExVol {

protected:

    Chain * chain;
    std::vector<BPStep*> BPS;
    arma::mat*  bp_pos;
    arma::cube* triads;

    int num_bp;
    int num_bps;
    double disc_len;

    bool closed_topology;

    bool check_crossings = true;

    int counter;
    int counter_reject;
    int counter_RP_reject;

    arma::mat  bp_pos_backup;
    arma::cube triads_backup;

    // EV_beads
    double EV_dist;            // actual exclusion radius
    double eff_size_EV_bead;  // the effective maximum size of an EV bead;
    int    num_bp_per_EV;
    int    upper_shift;
    int    num_EV;
    int    neighbour_skip;
    int    neighbour_skip_plus_one;
    int    neighbour_skip_boundary;
    int    neighbour_skip_boundary_plus_one;
    bool   check_beads;

    arma::ivec EV_beads;

    bool additional_boundcheck;
    std::vector<arma::ivec>  addboundpairs;

public:

    ExVol(Chain * ch, double EV_distance, bool check_crossings = true);
    ~ExVol();

    arma::mat*  get_bp_pos_backup();
    arma::cube* get_triads_backup();

    double get_EV_dist();
    int    get_bp_per_EV();


protected:

//    inline int bp2EV(int bp_id);
    inline int bp2EV_upper(int bp_id);
    inline int bp2EV_lower(int bp_id);

    arma::ivec interval_bp2EV(const arma::ivec& bp_interval);

    inline bool within_EV_typeA(int id);
    inline bool within_EV_typeB(int id);
    inline bool within_EV_typeC(int id);
    inline bool within_EV_typeD(int id);
    inline bool within_EV_typeE(int id);


public:
    bool check(const std::vector<arma::ivec>* moved);
    void set_backup_conf(arma::mat * backup_pos, arma::cube * backup_triads);

protected:
    void   cal_EV_intervals(    const std::vector<arma::ivec>* moved,
                                std::vector<arma::ivec>& EV_typeA,
                                std::vector<arma::ivec>& EV_typeB,
                                std::vector<arma::ivec>& EV_typeC,
                                std::vector<arma::ivec>& EV_typeD,
                                std::vector<arma::ivec>& EV_typeE);

    bool    check_intervals(    const std::vector<arma::ivec>& EV_typeA,
                                const std::vector<arma::ivec>& EV_typeB,
                                const std::vector<arma::ivec>& EV_typeC,
                                const std::vector<arma::ivec>& EV_typeD,
                                const std::vector<arma::ivec>& EV_typeE);

    bool    check_interval_singleMove(int A1, int A2, int B1, int B2);
    bool    check_interval_doubleMove(int A1, int A2, int B1, int B2);
    bool    check_within_interval(int A1, int A2);
    double  doubleMove(int id1, int id2);
    double  singleMove(int id1, int id2);


    bool    check_intervals_simpleoverlap(  const std::vector<arma::ivec>& EV_typeA,
                                            const std::vector<arma::ivec>& EV_typeB,
                                            const std::vector<arma::ivec>& EV_typeC,
                                            const std::vector<arma::ivec>& EV_typeD,
                                            const std::vector<arma::ivec>& EV_typeE);
    bool    check_interval_EVoverlap(int A1, int A2, int B1, int B2);
    double  beadoverlap(int id1, int id2);

    void    transfer_config(const std::vector<arma::ivec>* moved, arma::mat* bp_pos_from, arma::mat* bp_pos_to, arma::cube* triads_from, arma::cube* triads_to );

public:
    void    revert_to_backup();
    void    set_current_as_backup();

public:
    double rejection_rate();

/*
    Repulsion Plane Checks
*/
protected:
    bool repulsion_plane_active;
    int  first_last = REPULSIONPLANE_FIRST_AND_LAST;

public:
    void set_repulsion_plane(bool active=true);
    void set_repulsion_plane(bool plane_at_first, bool plane_at_last);

protected:
    bool RP_check_intervals(const std::vector<arma::ivec>& EV_typeB,
                            const std::vector<arma::ivec>& EV_typeC,
                            const std::vector<arma::ivec>& EV_typeD,
                            const std::vector<arma::ivec>& EV_typeE);
    bool RP_interval(int A, int B);

    bool RP_interval_first(int A, int B);
    bool RP_interval_last(int A, int B);

/*
    Line Closure
*/
protected:
    bool line_closure_active = false;
    bool line_closure_check_front = true;
    bool line_closure_check_back  = true;
    arma::colvec line_closure_direction =  {0,0,1};

public:
    void set_line_closure(bool active=true);
    void set_line_closure(bool check_front, bool check_back);

protected:
//    bool LC_check_intervals(const std::vector<arma::ivec>& EV_typeB,
//                            const std::vector<arma::ivec>& EV_typeC,
//                            const std::vector<arma::ivec>& EV_typeD,
//                            const std::vector<arma::ivec>& EV_typeE);
    bool LC_interval(int A, int B);


// Debug Methods
protected:
    bool check_moved_interval_order(const std::vector<arma::ivec>* moved);
public:
    bool check_overlap();


};


#endif
