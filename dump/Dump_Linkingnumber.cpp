#include "Dump_Linkingnumber.h"

Dump_Linkingnumber::Dump_Linkingnumber(Chain * ch, int N_dump, const std::string& filename, const std::string& options, bool append)
: Dump(ch,N_dump,filename,append)
{
    if (!app) {
        std::ofstream ofstr;
        ofstr.open(fn, std::ofstream::out | std::ofstream::trunc);
        ofstr.close();
    }

    iID = 0;
    fID = num_bp-1;

    if (options=="exact") {
        option=0;
    }
    if (options=="fuller") {
        option=1;
    }
    if (options=="both") {
        option=2;
    }
}
Dump_Linkingnumber::Dump_Linkingnumber(Chain * ch, int N_dump, const std::string& filename, int start_id, int end_id, const std::string& options, bool append)
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


    option = -1;
    if (options=="exact") {
        option=0;
    }
    if (options=="fuller") {
        option=1;
    }
    if (options=="both") {
        option=2;
    }
    if (options=="compare") {
        option=3;
    }
    if (options=="ebfu" || options=="fueb") {
        option=4;
    }

    if (option == -1) {
        std::cout << "Invalid option for Dump_Linkingnumber!" << std::endl;
        std::exit(0);
    }
}
Dump_Linkingnumber::~Dump_Linkingnumber() {}

void Dump_Linkingnumber::prod_dump() {

    double Tw = chain->cal_twist(iID,fID);
    if (option == 0) {
        arma::mat subpos = *chain->get_bp_pos();
        subpos           = subpos(arma::span(0,2),arma::span(iID,fID));
        double Wr        = chain->cal_langowski_writhe_1a(subpos,false);

//        double Wr = chain->cal_langowski_writhe_1a();
        std::ofstream ofstr;
        ofstr.open(fn, std::ofstream::out | std::ofstream::app);
        ofstr  << Tw << " " << Wr << "\n";
        ofstr.close();
    }
    if (option == 1) {
        arma::colvec force_dir = chain->get_force_dir();
        double Wr = chain->cal_fuller_writhe(iID,fID,force_dir);

//        double Wr = chain->cal_fuller_writhe(0,num_bp, force_dir);
        std::ofstream ofstr;
        ofstr.open(fn, std::ofstream::out | std::ofstream::app);
        ofstr  << Tw << " " << Wr << "\n";
        ofstr.close();
    }
    if (option == 2) {
        arma::colvec force_dir = chain->get_force_dir();
        double Wr_full = chain->cal_fuller_writhe(iID,fID, force_dir);

        arma::mat subpos = *chain->get_bp_pos();
        subpos           = subpos(arma::span(0,2),arma::span(iID,fID));
        double Wr_lang   = chain->cal_langowski_writhe_1a(subpos,false);

        double new_writhe = cal_fullerwrithe(subpos,force_dir);
        double closure    = cal_closure_writhe(subpos,force_dir);

        double Wr_lang_1a = cal_langowskiwrithe(subpos);
        double Wr_lang_1b = cal_langowskiwrithe_1b(subpos);

        std::ofstream ofstr;
        ofstr.open(fn, std::ofstream::out | std::ofstream::app);
//        ofstr  << Tw << " " << Wr_lang << " " << Wr_full << " " << Wr_lang_1a << " " << Wr_lang_1b << "\n";
        ofstr  << Tw << Wr_full << " " << Wr_lang_1b << "\n";
        ofstr.close();
    }

    if (option == 3) {
        arma::colvec force_dir = chain->get_force_dir();
//        double Wr_full = chain->cal_fuller_writhe(iID,fID, force_dir);



        arma::mat subpos = *chain->get_bp_pos();
//        subpos           = subpos(arma::span(0,2),arma::span(iID,fID));


//        double new_writhe = cal_fullerwrithe(subpos,force_dir);
//        double closure    = cal_closure_writhe(subpos,force_dir);


        int N = subpos.n_cols;
        int add_segs = 3;
        arma::mat extpos = arma::zeros(3,N+2*add_segs);

//        extpos.col(0) = subpos.col(0) - (subpos.col(1) - subpos.col(0));
//        extpos.col(N+1) = subpos.col(N-1) + subpos.col(N-1) - subpos.col(N-2);

        for (int i=add_segs;i<N+add_segs;i++) {
            extpos.col(i) = subpos.col(i-add_segs);
        }

        for (int i=0;i<add_segs;i++) {
            extpos.col(add_segs-1-i) = extpos.col(add_segs-i);
            extpos.col(add_segs-1-i)(2) -= disc_len;
        }
        for (int i=0;i<add_segs;i++) {
            extpos.col(N+add_segs+i) = extpos.col(N+add_segs-1+i);
            extpos.col(N+add_segs+i)(2) += disc_len;
        }


//        for (int i=add_segs-1;i>=0;i--) {
//            extpos.col(i) = extpos.col(i+1) + (extpos.col(i+1)-extpos.col(i+2));
//        }
//        for (int i=N+add_segs;i<N+2*add_segs;i++) {
//            extpos.col(i) = extpos.col(i-1) + (extpos.col(i-1)-extpos.col(i-2));
//        }


//        double Wr_lang_1a = cal_langowskiwrithe(extpos);
        double Wr_lang_1b = cal_langowskiwrithe_1b(extpos);

        double Wr_full = cal_fullerwrithe(extpos,force_dir);
//        double Wr_lang  = chain->cal_langowski_writhe_1a(extpos,false);

        double Lk_full = Tw + Wr_full;
//        double Lk_lang = Tw + Wr_lang;
        double Lk_bead = chain->get_dLK();
//        double Lk_lang_1a = Tw + Wr_lang_1a;
        double Lk_lang_1b = Tw + Wr_lang_1b;

        std::ofstream ofstr;
        ofstr.open(fn, std::ofstream::out | std::ofstream::app);
        ofstr  << Tw << " " << Lk_full << " " << Lk_lang_1b << " " << Lk_bead << "\n";
        ofstr.close();
    }

    if (option == 4) {
        arma::colvec force_dir = chain->get_force_dir();
        arma::mat subpos = *chain->get_bp_pos();
//        subpos           = subpos(arma::span(0,2),arma::span(iID,fID));

        double Wr_full = cal_fullerwrithe(subpos,force_dir);
        double Lk_full = Tw + Wr_full;
        double Lk_bead = chain->get_dLK();

        std::ofstream ofstr;
        ofstr.open(fn, std::ofstream::out | std::ofstream::app);
        ofstr  << Tw << " " << Lk_full << " " << Lk_bead << "\n";
        ofstr.close();
    }


    chain->get_dLK();
}


