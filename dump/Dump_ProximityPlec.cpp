#include "Dump_ProximityPlec.h"

Dump_ProximityPlec::Dump_ProximityPlec(Chain * ch, int N_dump, const std::string& filename, double threshold_dist, int min_seg_dist, double frac_closest, bool append)
: Dump(ch,N_dump,filename,append),
threshold_dist(threshold_dist),
min_seg_dist(min_seg_dist),
frac_closest(frac_closest)
{
    if (!app) {
        std::ofstream ofstr;
        ofstr.open(fn, std::ofstream::out | std::ofstream::trunc);
        ofstr.close();
    }
    proxplecfind = new ProximityPlecFinder(threshold_dist,min_seg_dist,frac_closest);
}
Dump_ProximityPlec::~Dump_ProximityPlec() {}


void Dump_ProximityPlec::prod_dump() {
    double xp = proxplecfind->get_xp(pos);

    std::ofstream ofstr;
    ofstr.open(fn, std::ofstream::out | std::ofstream::app);
    ofstr << xp << "\n";
    ofstr.close();
}

void Dump_ProximityPlec::final_dump() {
}







