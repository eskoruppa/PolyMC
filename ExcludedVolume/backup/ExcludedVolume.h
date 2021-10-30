#ifndef __EXCLUDED_VOLUME__
#define __EXCLUDED_VOLUME__

#include "../Chain.h"
#include "../SO3Methods.h"
#include "../ExtraFuncs.h"

#include "EV_Beads.h"
#include "ProximityCluster.h"

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

#define EV_FULLCHECK_VIOLATION 10000

/*
    DEBUG FLAG
*/

//#define DEBUG_EXVOL



/*
Different types of moved intervals:

TypeA:
    indicator range:        0
    type of displacement:   Beads that did not move and thus don't require crosschecking with other
                            non-moved beads

TypeB:
    indicator range:        1
    type of displacement:   Beads that were collectively trianslated. Require only single move check with
                            TypeA and Type D.

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


class ExVol;

class ExVol {

protected:

    Chain * chain;
    vector<BPStep*> BPS;
    arma::mat*  bp_pos;
    arma::cube* triads;

    int num_bp;
    int num_bps;
    double disc_len;

    int counter;
    int counter_reject;
    int counter_RP_reject;

    arma::mat  bp_pos_backup;
    arma::cube triads_backup;

    // EV_beads
    double EV_dist;            // actual exclusion radius
    double eff_size_EV_bead;  // the effective maximum size of an EV bead;
    int    num_bp_per_EV;
    int    num_EV;
    int    bpid_first_EV;

    arma::ivec EV_beads;
//    vector<ExcludedVolumeBead*> EV_beads;

    // ProximityCluster
    int     bpid_first_PC;
    int     EVid_first_PC;

    int    num_PC;
    int    num_EV_per_PC;
    int    num_bp_per_PC;
    double PC_radius;

    vector<ProximityCluster*> PC_list;

    long long int steps_since_last_PC_init;


public:

    ExVol(Chain * ch, double EV_radius, int EV_per_PC);
    ~ExVol();

    arma::mat*  get_bp_pos_backup();
    arma::cube* get_triads_backup();

    double get_EV_dist();

protected:

    inline int bp2EV(int bp_id);
    inline int bp2EV_upper(int bp_id);
    inline int bp2EV_lower(int bp_id);
    inline int EV2PC(int EV_id);
    inline int EV2PC_upper(int EV_id);
    inline int EV2PC_lower(int EV_id);
    inline int bp2PC(int bp_id);
    inline int bp2PC_upper(int bp_id);
    inline int bp2PC_lower(int bp_id);

    arma::ivec interval_bp2EV(const arma::ivec& bp_interval);
    arma::ivec interval_bp2PC(const arma::ivec& bp_interval);
    arma::ivec interval_EV2PC(const arma::ivec& EV_interval);

    inline bool within_EV_typeA(int id);
    inline bool within_EV_typeB(int id);
    inline bool within_EV_typeC(int id);
    inline bool within_EV_typeD(int id);
    inline bool within_EV_typeE(int id);


public:
    bool check(const vector<arma::ivec>* moved);
    void set_backup_conf(arma::mat * backup_pos, arma::cube * backup_triads);

protected:
    bool   update_PC(const vector<arma::ivec>* moved);
    void   cal_EV_intervals(    const vector<arma::ivec>* moved,
                                vector<arma::ivec>& EV_typeA,
                                vector<arma::ivec>& EV_typeB,
                                vector<arma::ivec>& EV_typeC,
                                vector<arma::ivec>& EV_typeD,
                                vector<arma::ivec>& EV_typeE);

    bool    check_intervals(    const vector<arma::ivec>& EV_typeA,
                                const vector<arma::ivec>& EV_typeB,
                                const vector<arma::ivec>& EV_typeC,
                                const vector<arma::ivec>& EV_typeD,
                                const vector<arma::ivec>& EV_typeE);

    bool    check_singleMove(int A1, int A2, int B1, int B2);
    bool    check_doubleMove(int A1, int A2, int B1, int B2);
    bool    check_within_interval(int A1, int A2);
    double  check_interval_boundary(int A, int B);
    double  doubleMove(int id1, int id2);
    double  singleMove(int id1, int id2);

    void    transfer_config(const vector<arma::ivec>* moved, arma::mat* bp_pos_from, arma::mat* bp_pos_to, arma::cube* triads_from, arma::cube* triads_to );

public:
    void    revert_to_backup();


public:
    double rejection_rate();


/*
    Repulsion Plane Checks
*/
protected:
    bool repulsion_plane_active;

public:
    void set_repulsion_plane(bool active=true);

protected:
    bool RP_check_intervals(const vector<arma::ivec>& EV_typeB,
                            const vector<arma::ivec>& EV_typeC,
                            const vector<arma::ivec>& EV_typeD,
                            const vector<arma::ivec>& EV_typeE);
    bool RP_interval(int A, int B);




// Testing
protected:
    bool check_moved_interval_order(const vector<arma::ivec>* moved);
    bool check_full();
public:
    bool check_overlap();






};

#endif
