#include "Dump_Energy.h"

/*
TODO:
    Dump stepwise energies. The boolian switch to activate this feature is included but doesn't do anything yet.
    This will require a new function in the BPS class that can extract the local energies without first selecting
    them for extraction by taking half of the long range coupling energies into account.
*/


Dump_Energy::Dump_Energy(Chain * ch, int N_dump, const std::string& filename, bool stepwise, bool append)
: Dump(ch,N_dump,filename,append),stepwise(stepwise)
{
    if (!app) {
        std::ofstream ofstr;
        ofstr.open(fn, std::ofstream::out | std::ofstream::trunc);
        ofstr.close();
    }
}
Dump_Energy::~Dump_Energy() {}


void Dump_Energy::prod_dump() {
    double energy = chain->extract_true_energy(); //*chain->get_T()/chain->get_T_ref();

    std::ofstream ofstr;
    ofstr.open(fn, std::ofstream::out | std::ofstream::app);
    ofstr << step << " " << energy;
//    if (stepwise) {
//        for (int i=0;i<num_bp;i++){
//            BPS[i]->
//        }
//    }
    ofstr << "\n";
    ofstr.close();
}
