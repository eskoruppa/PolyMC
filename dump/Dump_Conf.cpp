#include "Dump_Conf.h"

Dump_Conf::Dump_Conf(Chain * ch, int N_dump, const std::string& filename,bool append)
: Dump(ch,N_dump,filename,append)
{
    if (!app) {
        std::ofstream ofstr;
        ofstr.open(fn, std::ofstream::out | std::ofstream::trunc);
        ofstr.close();
    }

}
Dump_Conf::~Dump_Conf() {}

void Dump_Conf::prod_dump() {

    std::ofstream conffile;
    conffile.open(fn, std::ofstream::app);

    arma::vec translate = arma::sum(*chain->get_bp_pos(),1)*1.0/num_bp;

    arma::vec spos;
    for (int i=0;i<num_bp;i++) {
        spos = chain->get_bp_pos()->col(i)-translate;
        conffile << i+1 << " " << spos(0) << " " << spos(1) << " " << spos(2) << " \n";
    }

    conffile.close();
}

void Dump_Conf::final_dump() {
//    prod_dump();
}



