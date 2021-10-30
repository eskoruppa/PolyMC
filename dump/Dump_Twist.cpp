#include "Dump_Twist.h"

Dump_Twist::Dump_Twist(Chain * ch, int N_dump, const std::string& filename, bool append)
: Dump(ch,N_dump,filename,append)
{
    if (!app) {
        std::ofstream ofstr;
        ofstr.open(fn, std::ofstream::out | std::ofstream::trunc);
        ofstr.close();
    }
    iID = 0;
    fID = num_bp-1;
    init_store();
}
Dump_Twist::Dump_Twist(Chain * ch, int N_dump, const std::string& filename, int start_id, int end_id, bool append)
: Dump(ch,N_dump,filename,append)
{
    if (!app) {
        std::ofstream ofstr;
        ofstr.open(fn, std::ofstream::out | std::ofstream::trunc);
        ofstr.close();
    }
    iID = start_id;
    if (iID < 0) iID = 0;
    fID = end_id;
    if (fID > num_bp-2 || fID <= iID) fID = num_bp-2;
    init_store();
}

Dump_Twist::~Dump_Twist() {}

void Dump_Twist::init_store() {
    store_tw      = arma::zeros(DTW_STORE_TW);
    store_counter = 0;
}

void Dump_Twist::write2file() {
    std::ofstream ofstr;
    ofstr.open(fn, std::ofstream::out | std::ofstream::app);
    for (unsigned i=0;i<store_counter;i++) {
        ofstr << store_tw(i) << "\n";
    }
    ofstr.close();
    store_counter = 0;
}

void Dump_Twist::prod_dump() {
    store_tw(store_counter) = chain->cal_twist(iID,fID);
    store_counter++;
    if (store_counter == DTW_STORE_TW) {
        write2file();
    }
}

void Dump_Twist::final_dump() {
    write2file();
}

