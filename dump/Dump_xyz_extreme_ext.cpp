#include "Dump_xyz_extreme_ext.h"

Dump_xyz_extreme_ext::Dump_xyz_extreme_ext(Chain * ch, int N_dump, const std::string& filename, int min_datapoints, double num_sigmas,bool append, const std::string& center, const std::string& representation)
: Dump_xyz(ch,N_dump,filename,append),min_datapoints(min_datapoints),num_sigmas(num_sigmas), counter(0), center_conf(center), rep(representation)
{
    fn       = filename+".xyz";
    fn_stats = filename+".stats";

    if (!app) {
        std::ofstream ofstr;
        ofstr.open(fn, std::ofstream::out | std::ofstream::trunc);
        ofstr.close();
        ofstr.open(fn_stats, std::ofstream::out | std::ofstream::trunc);
        ofstr.close();
    }
    sum_z    = 0;
    sum_z_sq = 0;
    sum_num  = 0;
    cooldowncount = 0;

    num_snapshots = 0;
}
Dump_xyz_extreme_ext::~Dump_xyz_extreme_ext() {}

void Dump_xyz_extreme_ext::prod_dump() {
    double z = arma::dot(chain->get_force_dir(),chain->get_bp_pos()->col(num_bp-1)-chain->get_bp_pos()->col(0));
    sum_z    += z;
    sum_z_sq += z*z;
    sum_num++;
    if (cooldowncount > 0) {
        cooldowncount--;
    }
    if (sum_num > min_datapoints) {
        double mean  = sum_z/sum_num;
        double stdev = std::sqrt(sum_z_sq/sum_num - mean*mean);

        if ((z > mean + num_sigmas*stdev || z < mean - num_sigmas*stdev) && cooldowncount == 0) {
            std::cout << "range: [" <<  mean - num_sigmas*stdev << "," << mean + num_sigmas*stdev << "] - " << z << std::endl;
            xyz2file();
            cooldowncount = cooldownsteps;

            num_snapshots++;
            // dump 2 stats file
            double Tw = chain->cal_twist(0,num_bp);
            double Lk = chain->get_dLK();
            double Wr = Lk - Tw;

            std::ofstream ofstr;
            ofstr.open(fn_stats, std::ofstream::out | std::ofstream::app);
            ofstr << num_snapshots << " " << z << " " << mean << " " << sum_z_sq/sum_num << " " << Lk << " " << Tw << " " << Wr << std::endl;
            ofstr.close();
        }
    }

}

void Dump_xyz_extreme_ext::final_dump() {
//    prod_dump();
}



