#ifndef __MCS_BRANCHWINDER_H__
#define __MCS_BRANCHWINDER_H__
#include "MCStep.h"

#define MCSBW_HINGEFIND_RETRIES 10
#define MCSBW_MIN_HINGE_DIST    5

/*
    This variable controles whether B is always the closest bead to A or whether it
    is chosen randomly within the given distance constraint.
    TODO: Consider to remove the distance constraint
*/
#define MCSBW_CHOOSE_CLOSEST 1

/*
Perhaps exclude the possibility to rotate te first triad.

*/

class MCS_BranchWinder;

class MCS_BranchWinder: public MCStep
{
protected:
    bool    closed;
    int     n_beads_in_seg;
    int     min_contour_bead_dist;
    double  max_dist;

    double sigma;

    std::uniform_int_distribution<> genhingeA{0, 1};
    std::uniform_int_distribution<> selecthingeB{0, 1};

public:
    MCS_BranchWinder( Chain * ch, double seg_size, double max_dist, double min_contour_dist);
    ~MCS_BranchWinder();

    void update_settings();

// Monte Carlo Step (call main functionality)
public:
    bool MC_move();


protected:
    int     get_partner_random(int idA);
    int     get_partner_closest(int idA);
    int     get_partner_closest_test(int idA);
    bool    direction_is_forward(int idA,int idB);


};

#endif
