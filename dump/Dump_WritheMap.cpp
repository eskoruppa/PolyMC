#include "Dump_WritheMap.h"

Dump_WritheMap::Dump_WritheMap(Chain * ch, int N_dump, const std::string& filename, double seg_size, bool append)
: Dump(ch,N_dump,filename,append), seg_size(seg_size)
{
    if (!app) {
        std::ofstream ofstr;
        ofstr.open(fn, std::ofstream::out | std::ofstream::trunc);
        ofstr.close();
    }
}
Dump_WritheMap::~Dump_WritheMap() {}

void Dump_WritheMap::prod_dump() {
    std::ofstream file;
    file.open(fn, std::ofstream::app);

    arma::mat writhemat;
    chain->langowski_writhe_elements(&writhemat,seg_size);
    int elem = writhemat.n_cols;

    for (int i=0;i<elem;i++) {
        for (int j=0;j<elem-1;j++) {
            file << writhemat(i,j) << ",";
        }
        file << writhemat(i,elem-1) << "\n";
    }
    file.close();
}

void Dump_WritheMap::final_dump() {
//    prod_dump();
}



