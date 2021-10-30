#include "Dump_State.h"

Dump_State::Dump_State(Chain * ch, std::string mode, int N_dump, const std::string& filename,bool dump_triads,bool dump_Omegas,double EVrad,bool append)
: Dump(ch,N_dump,filename,append), counter(0), dump_triads(dump_triads), dump_Omegas(dump_Omegas), EV_included(false)
{
    if (!app || !fexists(fn)) {
        std::ofstream newfile;
        newfile.open(fn, std::ofstream::out | std::ofstream::trunc);
        newfile << "################################" << std::endl;
        newfile << "############ CONFIG ############" << std::endl;
        newfile << "pos     : " << dump_pos << std::endl;
        newfile << "triads  : " << dump_triads << std::endl;
        newfile << "Omegas  : " << dump_Omegas << std::endl;
        newfile << "Segments: " << chain->get_num_bp() << std::endl;
        newfile << "disc_len: " << chain->get_disc_len() << std::endl;
        newfile << "Ia_range: " << chain->get_interaction_range() << std::endl;
        newfile << "T       : " << chain->get_T() << std::endl;
        newfile << "closed  : " << chain->topology_closed() << std::endl;
        newfile << "Lk_fixed: " << chain->link_fixed() << std::endl;
        if (EVrad>0) {
            newfile << "EVactive: " << 1 << std::endl;
            newfile << "EVradius: " << EVrad << std::endl;
        }
        else {
            newfile << "EVactive: " << 0 << std::endl;
            newfile << "EVradius: " << 0 << std::endl;
        }
        newfile << "delta_LK: " << chain->get_dLK() << std::endl;
        newfile << "force   : " << chain->get_force() << std::endl;
        newfile << "################################" << std::endl;
        newfile.close();
    }
}


Dump_State::Dump_State(Chain * ch, std::string mode, int N_dump, const std::string& filename,bool dump_triads,bool dump_Omegas,bool append)
: Dump(ch,N_dump,filename,append), counter(0), dump_triads(dump_triads), dump_Omegas(dump_Omegas), EV_included(false)
{
    if (!app || !fexists(fn)) {
        std::ofstream newfile;
        newfile.open(fn, std::ofstream::out | std::ofstream::trunc);
        newfile << "################################" << std::endl;
        newfile << "############ CONFIG ############" << std::endl;
        newfile << "pos     : " << dump_pos << std::endl;
        newfile << "triads  : " << dump_triads << std::endl;
        newfile << "Omegas  : " << dump_Omegas << std::endl;
        newfile << "Mode    : " << mode << std::endl;
        newfile << "Segments: " << chain->get_num_bp() << std::endl;
        newfile << "disc_len: " << chain->get_disc_len() << std::endl;
        newfile << "Ia_range: " << chain->get_interaction_range() << std::endl;
        newfile << "T       : " << chain->get_T() << std::endl;
        newfile << "closed  : " << chain->topology_closed() << std::endl;
//        newfile << "Lk_fixed: " << chain->link_fixed() << std::endl;
        newfile << "EVactive: " << 0 << std::endl;
        newfile << "EVradius: " << 0 << std::endl;
        newfile << "delta_LK: " << chain->get_dLK() << std::endl;
        newfile << "force   : " << chain->get_force() << std::endl;
        newfile << "################################" << std::endl;
        newfile.close();
    }
}

Dump_State::Dump_State(Chain * ch, std::string mode, ExVol * EV, int N_dump, const std::string& filename,bool dump_triads,bool dump_Omegas,bool append)
: Dump(ch,N_dump,filename,append), EV(EV), EV_included(true), counter(0), dump_triads(dump_triads), dump_Omegas(dump_Omegas)
{
    if (!app || !fexists(fn)) {
        std::ofstream newfile;
        newfile.open(fn, std::ofstream::out | std::ofstream::trunc);
        newfile << "################################" << std::endl;
        newfile << "############ CONFIG ############" << std::endl;
        newfile << "pos     : " << dump_pos << std::endl;
        newfile << "triads  : " << dump_triads << std::endl;
        newfile << "Omegas  : " << dump_Omegas << std::endl;
        newfile << "Mode    : " << mode << std::endl;
        newfile << "Segments: " << chain->get_num_bp() << std::endl;
        newfile << "disc_len: " << chain->get_disc_len() << std::endl;
        newfile << "Ia_range: " << chain->get_interaction_range() << std::endl;
        newfile << "T       : " << chain->get_T() << std::endl;
        newfile << "closed  : " << chain->topology_closed() << std::endl;
//        newfile << "Lk_fixed: " << chain->link_fixed() << std::endl;
        newfile << "EVactive: " << 1 << std::endl;
        newfile << "EVradius: " << EV->get_EV_dist() << std::endl;
        newfile << "delta_LK: " << chain->get_dLK() << std::endl;
        newfile << "force   : " << chain->get_force() << std::endl;
        newfile << "################################" << std::endl;
        newfile.close();
    }
}

Dump_State::~Dump_State() {}

void Dump_State::prod_dump() {
    std::ofstream statefile;
    statefile.open(fn, std::ofstream::app);
    counter++;
    statefile << "Snapshot " << counter <<  " \n";

    if (dump_pos) {
        arma::vec spos;
        for (int i=0;i<num_bp;i++) {
            spos = chain->get_bp_pos()->col(i);
            statefile << std::setprecision(16) << spos(0) << " " << spos(1) << " " << spos(2) << " \n";
        }
    }

    if (dump_triads) {
        arma::mat triad;
        for (int i=0;i<num_bp;i++) {
            triad = chain->get_triads()->slice(i);
            statefile << std::setprecision(16) << triad.col(0)(0) << " " << triad.col(0)(1) << " " << triad.col(0)(2) << " "  << triad.col(1)(0) << " " << triad.col(1)(1) << " " << triad.col(1)(2) << " "  << triad.col(2)(0) << " " << triad.col(2)(1) << " " << triad.col(2)(2) << " \n";
        }
    }

    if (dump_Omegas) {
        arma::vec Omega;
        for (int i=0;i<num_bps;i++) {
            Omega = *BPS[i]->get_Theta();
            statefile << std::setprecision(16) << Omega(0) << " " << Omega(1) << " " << Omega(2) << " \n";
        }
    }
    statefile.close();
//    for (int i=0;i<1;i++) {
//        std::cout << i << std::endl;
//        PlecFinder pf(chain);
//    }
}

void Dump_State::final_dump() {
//    prod_dump();
}
















































