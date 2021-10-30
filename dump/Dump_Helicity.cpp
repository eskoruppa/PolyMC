#include "Dump_Helicity.h"

Dump_Helicity::Dump_Helicity(Chain * ch, int N_dump, int N_dump2file, const std::string& filename, bool append)
: Dump(ch,N_dump,filename,append),
N_dump2file(N_dump2file),
dump2file_counter(0),
hel_sum(0)
{
    if (!app) {
        std::ofstream ofstr;
        ofstr.open(fn, std::ofstream::out | std::ofstream::trunc);
        ofstr.close();
    }
}
Dump_Helicity::~Dump_Helicity() {}

void Dump_Helicity::prod_dump() {

    double hel = 0;
    for (int i=0;i<num_bp-1;i++) {
        hel += arma::dot(arma::cross(triads->slice(i),triads->slice(i+1)),chain->get_force_dir());
    }
    hel = hel/(4*M_PI);
    hel_sum += hel;
    dump2file_counter++;

    if (dump2file_counter%N_dump2file==0) {
        std::ofstream ofstr;
        ofstr.open(fn, std::ofstream::out | std::ofstream::app);
        ofstr  << chain->get_dLK() << " " << hel_sum/dump2file_counter << std::endl;
        ofstr.close();

        hel_sum           = 0;
        dump2file_counter = 0;
    }
}


