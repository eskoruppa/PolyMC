#ifndef __DUMP_PLASMIDENDPOINTS_H__
#define __DUMP_PLASMIDENDPOINTS_H__
#include "Dump.h"

// Length of the partial segments in nm
#define PP_LOOP_LENGTH      340 //240//340
// Smooth distance in nm
#define PP_SMOOTHBOX        5 //50
// Minimum distance in nm between individual peaks in the spectrum of local writhe
#define PP_MIN_DETECT_DIST  80

// Peak detection theshold
#define PP_PEAK_THRESHOLD   1.0 //1.15
// Theshold writhe difference to distinguish peaks
#define PP_LOCAL_MIN_DIFF   0.4
// Agreement rate for crosschecks
#define PP_MIN_AGREEMENT    0.5  // 50%



class Dump_PlasmidEndpoints;

class Dump_PlasmidEndpoints: public Dump
{
protected:
    double density;
    int num_points;
    int point_dist;
    int loop_len;
    int half_length;
    int sample_step_dist;
    int sample_points;

    int  call_count     = 0;
    int  sample_count   = 0;
    bool sampling       = false;

    bool dump_xyz;
    bool dump_writhe;

    int  min_peak_dist;

//    arma::ivec sample_endpoints;


public:
    Dump_PlasmidEndpoints(Chain * ch, int N_dump, const std::string& filename, double density,int sample_step_dist, int sample_points,bool dump_writhe, bool print_xyz, bool append=false);
    ~Dump_PlasmidEndpoints();

    void prod_dump();
    void final_dump();

protected:

    arma::colvec smooth(arma::colvec& write, int box_size);
    int          find_endpoint_agreement(arma::ivec& num_endpoints);

    std::vector<arma::colvec> find_endpoints(arma::colvec & smoothed);
    std::vector<arma::colvec> find_endpoints_old(arma::colvec & smoothed);
};

#endif
