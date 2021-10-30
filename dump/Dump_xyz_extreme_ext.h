#ifndef __DUMP_XYZ_EXTREME_EXT_H__
#define __DUMP_XYZ_EXTREME_EXT_H__
#include "Dump.h"
#include "Dump_xyz.h"


class Dump_xyz_extreme_ext;

class Dump_xyz_extreme_ext: public Dump_xyz
{
protected:
    int counter;
    std::string center_conf;
    std::string rep;
    int min_datapoints;
    double num_sigmas;

    std::string fn_stats; // filename


    double sum_z;
    double sum_z_sq;
    int    sum_num;

    int    cooldowncount=0;
    int    cooldownsteps=100000;

    int    num_snapshots;


public:
    Dump_xyz_extreme_ext(Chain * ch, int N_dump, const std::string& filename, int min_datapoints, double num_sigmas, bool append=false, const std::string& center="COM", const std::string& representation="simple");
    ~Dump_xyz_extreme_ext();

    void prod_dump();
    void final_dump();



};




#endif
