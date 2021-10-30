#include "Dump_DistanceMap.h"

Dump_DistMap::Dump_DistMap(Chain * ch, int N_dump, const std::string& filename, double seg_size, bool append)
: Dump(ch,N_dump,filename,append)
{
    density = chain->get_disc_len()/seg_size;
    if (density>1) density=1;

    if (!app) {
        std::ofstream ofstr;
        ofstr.open(fn, std::ofstream::out | std::ofstream::trunc);
        ofstr.close();

        ofstr.open("dump/branch_test/conf.xyz", std::ofstream::out | std::ofstream::trunc);
        ofstr.close();
    }
}
Dump_DistMap::~Dump_DistMap() {}

void Dump_DistMap::prod_dump() {

    std::ofstream file;
    file.open(fn, std::ofstream::app);
    arma::mat dists = chain->distmap(density);
    int elem = dists.n_cols;

    for (int i=0;i<elem;i++) {
        for (int j=0;j<elem-1;j++) {
            file << dists(i,j) << ",";
        }
        file << dists(i,elem-1) << "\n";
    }
    file.close();


    std::ofstream xyzfile;
    xyzfile.open("dump/branch_test/conf.xyz", std::ofstream::app);
    xyzfile << num_bp << " \n";
    xyzfile << "Atoms. Timestep: " << step << " \n";
    arma::vec spos;
    for (int i=0;i<num_bp;i++) {
        spos = chain->get_bp_pos()->col(i);
        xyzfile << "A " <<  spos(0) << " " << spos(1) << " " << spos(2) << " \n";
    }

}

void Dump_DistMap::final_dump() {
//    prod_dump();
}



