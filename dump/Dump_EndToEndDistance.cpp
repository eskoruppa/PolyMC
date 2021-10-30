#include "Dump_EndToEndDistance.h"

Dump_EndToEndDistance::Dump_EndToEndDistance(Chain * ch, int N_dump, const std::string& filename, bool append)
: Dump(ch,N_dump,filename,append)
{
    if (!app) {
        std::ofstream ofstr;
        ofstr.open(fn, std::ofstream::out | std::ofstream::trunc);
        ofstr.close();
    }
}
Dump_EndToEndDistance::~Dump_EndToEndDistance() {}

void Dump_EndToEndDistance::prod_dump() {
    std::ofstream ofstr;
    ofstr.open(fn, std::ofstream::out | std::ofstream::app);

    ofstr << arma::norm(chain->get_bp_pos()->col(num_bp-1)-chain->get_bp_pos()->col(0)) << "\n";
    ofstr.close();

}

void Dump_EndToEndDistance::final_dump() {
}

