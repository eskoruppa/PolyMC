#include "Dump_Restart.h"

Dump_Restart::Dump_Restart(Chain * ch, int N_dump, const std::string& filename, bool append, bool last_snapshot)
: Dump(ch,N_dump,filename,append), counter(0), last(last_snapshot)
{
    if (!app && !last) {
        std::ofstream ofstr;
        ofstr.open(fn, std::ofstream::out | std::ofstream::trunc);
        ofstr.close();
    }
}


Dump_Restart::~Dump_Restart() {}

void Dump_Restart::prod_dump() {
    std::ofstream restartfile;
    restartfile.open(fn, std::ofstream::app);
    counter++;
    restartfile << "###############################################################" << "\n";
    restartfile << "snapshot: " << counter << " \n";
    restartfile << "type:     " << chain->get_config_type() << " \n";
    restartfile << "num_bp:   " << chain->get_num_bp() << " \n";
    restartfile << "sequence: " << chain->get_sequence() << " \n";
    restartfile << "dLK:      " << chain->get_dLK() << " \n";

    arma::vec spos;
    for (int i=0;i<num_bp;i++) {
        spos = chain->get_bp_pos()->col(i);
        restartfile << std::setprecision(14) << spos(0) << " " << spos(1) << " " << spos(2) << " \n";
    }

    arma::mat triad;
    for (int i=0;i<num_bp;i++) {
        triad = chain->get_triads()->slice(i);
        restartfile << std::setprecision(14) << triad.col(0)(0) << " " << triad.col(0)(1) << " " << triad.col(0)(2) << " "  << triad.col(1)(0) << " " << triad.col(1)(1) << " " << triad.col(1)(2) << " "  << triad.col(2)(0) << " " << triad.col(2)(1) << " " << triad.col(2)(2) << " \n";
    }
}

void Dump_Restart::final_dump() {
    if (last) {
        std::string dupfn = fn;
        int dupcount = 1;
        while (fexists(dupfn)) {
            dupcount++;
            dupfn = fn + "_#" + std::to_string(dupcount);
        }
        fn = dupfn;
        std::ofstream ofstr;
        ofstr.open(fn, std::ofstream::out | std::ofstream::trunc);
        ofstr.close();
        counter = -2;
        prod_dump();
    }
//    else {
//        prod_dump();
//    }
}

