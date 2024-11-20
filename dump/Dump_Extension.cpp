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
    multi_store = false;
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
    if (start_id < 0) {
        start_id = 0;
    }
    if (end_id < 0) {
        end_id = num_bp-1;
    }
    iID = start_id;
    if (iID < 0) iID = 0;
    fID = end_id;
    if (fID > num_bp-1 || fID <= iID) fID = num_bp-1;
    multi_store = false;
    init_store();
}

Dump_Extension::Dump_Extension(Chain * ch, int N_dump, const std::string& filename, int start_id, int end_id, int subdomain_size, bool append)
: Dump(ch,N_dump,filename,append)
{
    if (!app) {
        std::ofstream ofstr;
        ofstr.open(fn, std::ofstream::out | std::ofstream::trunc);
        ofstr.close();
    }
    if (start_id < 0) {
        start_id = 0;
    }
    if (end_id < 0) {
        end_id = num_bp-1;
    }
    iID = start_id;
    if (iID < 0) iID = 0;
    fID = end_id;
    if (fID > num_bp-1 || fID <= iID) fID = num_bp-1;
    multi_store = true;
    multi_num = std::ceil((fID-iID)/subdomain_size);
    multi_ids = arma::zeros(multi_num+1);
    // define ids
    for (unsigned i=0;i<multi_num;i++){
        multi_ids(i) = iID + i*subdomain_size;
    }
    multi_ids(multi_num) = fID;
    init_store();
}


Dump_Extension::~Dump_Extension() {}

void Dump_Extension::init_store() {
    if (multi_store){
        multi_store_exts = arma::zeros(DEXT_STORE_EXTS,multi_num);
    }
    else {
        store_exts    = arma::zeros(DEXT_STORE_EXTS);
    }
    store_counter = 0;
}

void Dump_Extension::write2file() {
    std::ofstream ofstr;
    ofstr.open(fn, std::ofstream::out | std::ofstream::app);
    if (multi_store) {
        for (unsigned i=0;i<store_counter;i++) {
            for (unsigned j=0;j<multi_num-1;j++) {
                // std::cout << i << " " << j << std::endl;
                // std::cout << multi_store_exts[i,j] << std::endl; 
                // std::cout << "printed" << std::endl;
                ofstr << multi_store_exts(i,j) << " ";
            }
            ofstr << multi_store_exts(i,multi_num-1) << "\n";
        }
    }
    else {
        for (unsigned i=0;i<store_counter;i++) {
            ofstr << store_exts(i) << "\n";
        }
    }
    ofstr.close();
    store_counter = 0;
}

void Dump_Extension::prod_dump() {
    if (multi_store) {
        for (unsigned i=0;i<multi_num;i++){
            multi_store_exts(store_counter,i) = arma::dot(chain->get_force_dir(),chain->get_bp_pos()->col(multi_ids(i+1))-chain->get_bp_pos()->col(multi_ids(i)));
        }
    }
    else {
        store_exts(store_counter) = arma::dot(chain->get_force_dir(),chain->get_bp_pos()->col(fID)-chain->get_bp_pos()->col(iID));
    }
    store_counter++;
    if (store_counter == DEXT_STORE_EXTS) {
        write2file();
    }
}

void Dump_Extension::final_dump() {
    write2file();
}

