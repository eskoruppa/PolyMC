#include "Dump_BubbleStats.h"

Dump_BubbleStats::Dump_BubbleStats(Chain * ch, int N_dump, const std::string& filename, bool append)
: Dump(ch,N_dump,filename,append)
{
    if (!app) {
        std::ofstream ofstr;
        ofstr.open(fn, std::ofstream::out | std::ofstream::trunc);
        ofstr.close();
    }
    // iID = 0;
    // fID = num_bp-1;
    init_store();
}
// Dump_BubbleStats::Dump_BubbleStats(Chain * ch, int N_dump, const std::string& filename, int start_id, int end_id, bool append)
// : Dump(ch,N_dump,filename,append)
// {
//     if (!app) {
//         std::ofstream ofstr;
//         ofstr.open(fn, std::ofstream::out | std::ofstream::trunc);
//         ofstr.close();
//     }
//     iID = start_id;
//     if (iID < 0) iID = 0;
//     fID = end_id;
//     if (fID > num_bp-1 || fID <= iID) fID = num_bp-1;
//     init_store();
// }

Dump_BubbleStats::~Dump_BubbleStats() {}

void Dump_BubbleStats::init_store() {
    store_num_bubsegs = arma::zeros(DBUB_STORE_NUM);
    store_num_bubs    = arma::zeros(DBUB_STORE_NUM);
    store_counter = 0;
}

void Dump_BubbleStats::write2file() {
    std::ofstream ofstr;
    ofstr.open(fn, std::ofstream::out | std::ofstream::app);
    for (unsigned i=0;i<store_counter;i++) {
        ofstr << store_num_bubsegs(i) << " " << store_num_bubs(i) << "\n";
    }
    ofstr.close();
    store_counter = 0;
}

void Dump_BubbleStats::prod_dump() {

    int num_bubsegs = arma::sum((*chain->get_states()))- chain->get_states()->size();
    int num_bubs = 0;
    if (chain->get_states()->at(0) == 2 ) {
        num_bubs += 1;
    }
    for (int i=0;i<chain->get_states()->size()-1;i++) {
        if (chain->get_states()->at(i) == 1 && chain->get_states()->at(i+1) == 2) {
            num_bubs += 1;
        }
    }
    store_num_bubsegs(store_counter) = num_bubsegs;
    store_num_bubs(store_counter)    = num_bubs;
    store_counter++;
    if (store_counter == DBUB_STORE_NUM) {
        write2file();
    }
}

void Dump_BubbleStats::final_dump() {
    write2file();
}

