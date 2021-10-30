#include "Dump_Extension.h"

Dump_Extension::Dump_Extension(Chain * ch, int N_dump, const std::string& filename, bool append)
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
Dump_Extension::Dump_Extension(Chain * ch, int N_dump, const std::string& filename, int start_id, int end_id, bool append)
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
    if (fID > num_bp-1 || fID <= iID) fID = num_bp-1;
    init_store();
}

Dump_Extension::~Dump_Extension() {}

void Dump_Extension::init_store() {
    store_exts    = arma::zeros(DEXT_STORE_EXTS);
    store_counter = 0;
}

void Dump_Extension::write2file() {
    std::ofstream ofstr;
    ofstr.open(fn, std::ofstream::out | std::ofstream::app);
    for (unsigned i=0;i<store_counter;i++) {
        ofstr << store_exts(i) << "\n";
    }
    ofstr.close();
    store_counter = 0;
}

void Dump_Extension::prod_dump() {
    store_exts(store_counter) = arma::dot(chain->get_force_dir(),chain->get_bp_pos()->col(fID)-chain->get_bp_pos()->col(iID));
    store_counter++;
    if (store_counter == DEXT_STORE_EXTS) {
        write2file();
    }
}

void Dump_Extension::final_dump() {
    write2file();
}

