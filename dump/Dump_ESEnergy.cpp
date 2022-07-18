#include "Dump_ESEnergy.h"


Dump_ESEnergy::Dump_ESEnergy(Chain * ch, int N_dump, const std::string& filename, ElStat * es,bool kt_units,bool recal, bool append)
: Dump(ch,N_dump,filename,append), elstat(es),kt_units(kt_units),recal(recal)
{
    if (!app) {
        std::ofstream ofstr;
        ofstr.open(fn, std::ofstream::out | std::ofstream::trunc);
        ofstr.close();
    }
}
Dump_ESEnergy::~Dump_ESEnergy() {}


void Dump_ESEnergy::prod_dump() {

    double energy = elstat->get_current_energy(kt_units,recal);

    std::ofstream ofstr;
    ofstr.open(fn, std::ofstream::out | std::ofstream::app);
    ofstr << step << " " << energy;
    ofstr << "\n";
    ofstr.close();
}
