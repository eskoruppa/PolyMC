#include "Dump_xyz.h"

Dump_xyz::Dump_xyz(Chain * ch, int N_dump, const std::string& filename,bool append, const std::string& center, const std::string& representation)
: Dump(ch,N_dump,filename,append), counter(0), center_conf(center), rep(representation)
{
    if (!app) {
        std::ofstream ofstr;
        ofstr.open(fn, std::ofstream::out | std::ofstream::trunc);
        ofstr.close();
    }

}
Dump_xyz::~Dump_xyz() {}

void Dump_xyz::prod_dump() {
    xyz2file();
}

void Dump_xyz::final_dump() {
//    prod_dump();
}


void Dump_xyz::xyz2file() {
    std::ofstream xyzfile;
    xyzfile.open(fn, std::ofstream::app);

    arma::vec translate = arma::zeros(3);
    if (center_conf=="first" || center_conf=="First") translate = chain->get_bp_pos()->col(0);
    if (center_conf=="com"   || center_conf=="COM")   translate = arma::sum(*chain->get_bp_pos(),1)*1.0/num_bp;
    if (center_conf=="last"  || center_conf=="Last")  translate = chain->get_bp_pos()->col(num_bp-1);

    if (rep == "simple") {
        xyzfile << num_bp << " \n";
        xyzfile << "Atoms. Timestep: " << step << " \n";
        arma::vec spos;
        for (int i=0;i<num_bp;i++) {
            spos = chain->get_bp_pos()->col(i)-translate;
            xyzfile << "A " <<  spos(0) << " " << spos(1) << " " << spos(2) << " \n";
        }
    }

    if (rep == "fl") {
        xyzfile << num_bp << " \n";
        xyzfile << "Atoms. Timestep: " << step << " \n";
        arma::vec spos;
        spos = chain->get_bp_pos()->col(0)-translate;
        xyzfile << "N " <<  spos(0) << " " << spos(1) << " " << spos(2) << " \n";
        for (int i=1;i<num_bp-1;i++) {
            spos = chain->get_bp_pos()->col(i)-translate;
            xyzfile << "A " <<  spos(0) << " " << spos(1) << " " << spos(2) << " \n";
        }
        spos = chain->get_bp_pos()->col(num_bp-1)-translate;
        xyzfile << "C " <<  spos(0) << " " << spos(1) << " " << spos(2) << " \n";
    }

    if (rep == "helix" || rep == "Helix") {
        xyzfile << num_bp*2 << " \n";
        xyzfile << "Atoms. Timestep: " << step << " \n";
        arma::vec spos;
        arma::vec bbpos;
        double dist = 0.5;
        for (int i=0;i<num_bp;i++) {
            spos = chain->get_bp_pos()->col(i)-translate;
            //~ xyzfile << "A " <<  spos(0) << " " << spos(1) << " " << spos(2) << " \n";
            bbpos = spos + chain->get_triads()->slice(i).col(1)*dist;
            xyzfile << "C " <<  bbpos(0) << " " << bbpos(1) << " " << bbpos(2) << " \n";
            bbpos = spos - chain->get_triads()->slice(i).col(1)*dist;
            xyzfile << "N " <<  bbpos(0) << " " << bbpos(1) << " " << bbpos(2) << " \n";
        }
    }
    if (rep == "dna" || rep == "DNA" ) {
        xyzfile << num_bp*4 << " \n";
        xyzfile << "Atoms. Timestep: " << step << " \n";
        arma::vec spos;
        arma::vec bbpos;
        double dist     = 0.5;
        double inc      = 0.15;
        double basedist = 0.1;
        for (int i=0;i<num_bp;i++) {
            spos = chain->get_bp_pos()->col(i)-translate;
            //~ xyzfile << "A " <<  spos(0) << " " << spos(1) << " " << spos(2) << " \n";
            bbpos = spos + chain->get_triads()->slice(i).col(1)*dist-chain->get_triads()->slice(i).col(0)*inc;
            xyzfile << "C " <<  bbpos(0) << " " << bbpos(1) << " " << bbpos(2) << " \n";
            bbpos = spos - chain->get_triads()->slice(i).col(1)*dist-chain->get_triads()->slice(i).col(0)*inc;
            xyzfile << "N " <<  bbpos(0) << " " << bbpos(1) << " " << bbpos(2) << " \n";
            bbpos = spos +   chain->get_triads()->slice(i).col(1)*basedist;
            xyzfile << "A " <<  bbpos(0) << " " << bbpos(1) << " " << bbpos(2) << " \n";
            bbpos = spos - chain->get_triads()->slice(i).col(1)*basedist;
            xyzfile << "A " <<  bbpos(0) << " " << bbpos(1) << " " << bbpos(2) << " \n";
        }
    }
    if (rep == "simpledna" || rep == "simpleDNA" ) {
        xyzfile << num_bp*2 << " \n";
        xyzfile << "Atoms. Timestep: " << step << " \n";
        arma::vec spos;
        arma::vec bbpos;
        double dist = 0.5;
        double inc  = 0.15;
        for (int i=0;i<num_bp;i++) {
            spos = chain->get_bp_pos()->col(i)-translate;
            //~ xyzfile << "A " <<  spos(0) << " " << spos(1) << " " << spos(2) << " \n";
            bbpos = spos + chain->get_triads()->slice(i).col(1)*dist-chain->get_triads()->slice(i).col(0)*inc;
            xyzfile << "C " <<  bbpos(0) << " " << bbpos(1) << " " << bbpos(2) << " \n";
            bbpos = spos - chain->get_triads()->slice(i).col(1)*dist-chain->get_triads()->slice(i).col(0)*inc;
            xyzfile << "N " <<  bbpos(0) << " " << bbpos(1) << " " << bbpos(2) << " \n";
        }
    }
    if (rep == "pdb" || rep == "PDB" ) {
        arma::vec line;
        arma::vec pos;
        arma::mat triad;
        for (int i=0;i<num_bp;i++) {
            triad = chain->get_triads()->slice(i);
            triad = triad.t();
            for (int j=0;j<3;j++) {
                line = triad.col(j);
                xyzfile <<  line(0) << " " << line(1) << " " << line(2) << "\n";
            }
        }
        for (int i=0;i<num_bp;i++) {
            pos   = chain->get_bp_pos()->col(i)-translate;
            xyzfile <<  pos(0) << " " << pos(1) << " " << pos(2) << "\n";
        }
    }
    if (rep == "triad" || rep == "Triad" ) {
        xyzfile << num_bp*2 << " \n";
        xyzfile << "Atoms. Timestep: " << step << " \n";
        arma::vec spos;
        arma::vec bbpos;
        double dist = 0.5;
        for (int i=0;i<num_bp;i++) {
            spos = chain->get_bp_pos()->col(i)-translate;
            xyzfile << "A " <<  spos(0) << " " << spos(1) << " " << spos(2) << " \n";
            bbpos = spos - chain->get_triads()->slice(i).col(1)*dist;
            xyzfile << "N " <<  bbpos(0) << " " << bbpos(1) << " " << bbpos(2) << " \n";
        }
    }
    if (rep == "triadf" || rep == "Triad_Full" ) {
        xyzfile << num_bp*3 << " \n";
        xyzfile << "Atoms. Timestep: " << step << " \n";
        arma::vec spos;
        arma::vec bbpos;
        double dist = 0.5;
        for (int i=0;i<num_bp;i++) {
            spos = chain->get_bp_pos()->col(i)-translate;
            xyzfile << "A " <<  spos(0) << " " << spos(1) << " " << spos(2) << " \n";
            bbpos = spos + chain->get_triads()->slice(i).col(0)*dist;
            xyzfile << "C " <<  bbpos(0) << " " << bbpos(1) << " " << bbpos(2) << " \n";
            bbpos = spos - chain->get_triads()->slice(i).col(1)*dist;
            xyzfile << "N " <<  bbpos(0) << " " << bbpos(1) << " " << bbpos(2) << " \n";
        }
    }
    if (rep == "EV" ) {
        int num_EV = 0;
        for (int i=2;i<num_bp;i=i+5) {
            num_EV++;
        }
        xyzfile << num_EV << " \n";
        xyzfile << "Atoms. Timestep: " << step << " \n";
        arma::vec spos;
        arma::vec bbpos;

        spos = chain->get_bp_pos()->col(2)-translate;
        xyzfile << "N " <<  spos(0) << " " << spos(1) << " " << spos(2) << " \n";

        int c = 0;
        for (int i=2+5;i<num_bp;i=i+5) {
            spos = chain->get_bp_pos()->col(i)-translate;

            if (c%3==0) {
                xyzfile << "A " <<  spos(0) << " " << spos(1) << " " << spos(2) << " \n";
            }
            if (c%3==1) {
                xyzfile << "C " <<  spos(0) << " " << spos(1) << " " << spos(2) << " \n";
            }
            if (c%3==2) {
                xyzfile << "N " <<  spos(0) << " " << spos(1) << " " << spos(2) << " \n";
            }
            c++;
        }
    }
    xyzfile.close();
}

