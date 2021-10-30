#ifndef __DUMP_PROXIMITYPLEC_H__
#define __DUMP_PROXIMITYPLEC_H__
#include "Dump.h"
#include "../PlecFinder/ProximityPlecFinder.h"

////////////////////////////////////////////////////////////////////////////////
/*
    Dump the extension of Writhe of a Tweezer Simulations
    Thus far only works in the constant linking number ensemble
*/

class Dump_ProximityPlec;
class Dump_ProximityPlec: public Dump
{
protected:
    double threshold_dist;
    int    min_seg_dist;
    double frac_closest;

    ProximityPlecFinder * proxplecfind;

public:
    Dump_ProximityPlec(Chain * ch, int N_dump, const std::string& filename, double threshold_dist, int min_seg_dist,double frac_closest, bool append=true);
    ~Dump_ProximityPlec();

    void prod_dump();
    void final_dump();
};

#endif
