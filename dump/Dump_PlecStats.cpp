#include "Dump_PlecStats.h"

Dump_PlecStats::Dump_PlecStats(Chain * ch, int N_dump, const std::string& filename, bool append)
: Dump(ch,N_dump,filename,append)
{
    if (!app) {
        std::ofstream ofstr;
        ofstr.open(fn, std::ofstream::out | std::ofstream::trunc);
        ofstr.close();
    }

}
Dump_PlecStats::~Dump_PlecStats() {}


void Dump_PlecStats::prod_dump() {
    PlecFinder plecs(chain);
    std::vector<double> stats = plecs.getPlecStats();

    std::ofstream ofstr;
    ofstr.open(fn, std::ofstream::out | std::ofstream::app);

    ofstr << stats[0];
    for (unsigned i=1;i<stats.size();i++) {
        ofstr << " " << stats[i];
    }
    ofstr << "\n";
    ofstr.close();
}

void Dump_PlecStats::final_dump() {
}







