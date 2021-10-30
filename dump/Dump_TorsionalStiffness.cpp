#include "Dump_TorsionalStiffness.h"

Dump_EffTorsStiff::Dump_EffTorsStiff(Chain * ch, int N_dump, const std::string& filename, const arma::colvec& force_dir, const std::string& writhe_method, int start_id, int end_id, bool append)
: Dump(ch,N_dump,filename,append)
{
    dLK = 0, dLK_sq = 0, Wr = 0, Wr_sq = 0, Tw = 0, Tw_sq = 0;
    counter = 0;

    iID = start_id;
    if (iID < 0) iID = 0;
    fID = end_id;
    if (fID > num_bp-1 || fID <= iID) fID = num_bp-1;

    method  = writhe_method;
    method_id = ETS_UNDEFINED;
    if (method=="fuller")       method_id = ETS_FULLER;
    if (method=="langowski")    method_id = ETS_LANGOWSKI;
    if (method=="both")         method_id = ETS_BOTH;

    if (method_id == ETS_UNDEFINED) {
        std::cout << "Error: Dump_EffTorsStiff:Dump_EffTorsStiff - Invalid method specified!" << std::endl;
        std::exit(0);
    }

    dir = force_dir;
}
Dump_EffTorsStiff::~Dump_EffTorsStiff() {}


void Dump_EffTorsStiff::prod_dump() {
    double twist = chain->cal_twist(iID,fID);
    double writhe;
    double link;

    if (method_id== ETS_FULLER) {
        writhe = chain->cal_fuller_writhe(iID,fID,dir);
        Wr      += writhe;
        Wr_sq   += writhe*writhe;
    }
    if (method_id== ETS_LANGOWSKI) {
        arma::mat subpos = *chain->get_bp_pos();
        subpos = subpos(arma::span(0,2),arma::span(iID,fID));
        writhe = chain->cal_langowski_writhe_1a(subpos,false);
//        writhe = chain->cal_langowski_writhe_1a();
        Wr      += writhe;
        Wr_sq   += writhe*writhe;
    }
    if (method_id== ETS_BOTH) {
        arma::mat subpos = *chain->get_bp_pos();
        subpos = subpos(arma::span(0,2),arma::span(iID,fID));
        double lanwr = chain->cal_langowski_writhe_1a(subpos,false);
//        writhe = chain->cal_langowski_writhe_1a();
        Wr      += lanwr;
        Wr_sq   += lanwr*lanwr;

        writhe     = chain->cal_fuller_writhe(iID,fID,dir);
        Wr_ful    += writhe;
        Wr_sq_ful += writhe*writhe;
    }

    link     = twist+writhe;
    dLK     += link;
    dLK_sq  += link*link;
    Tw      += twist;
    Tw_sq   += twist*twist;
    counter++;
}


void Dump_EffTorsStiff::final_dump() {
    dLK     /= counter;
    dLK_sq  /= counter;
    Wr      /= counter;
    Wr_sq   /= counter;
    Tw      /= counter;
    Tw_sq   /= counter;

    if (method_id== ETS_BOTH) {
        Wr_ful    /= counter;
        Wr_sq_ful /= counter;
    }


    double L    = chain->get_disc_len()*(fID-iID-1);
    double Ceff = L/(4*M_PI*M_PI)/(dLK_sq-dLK*dLK);

	std::ofstream ofstr;
	if (app) ofstr.open(fn, std::ofstream::out | std::ofstream::app);
	else     ofstr.open(fn, std::ofstream::out | std::ofstream::trunc);

    if (method_id== ETS_BOTH) {
        ofstr << chain->get_force() << " " << Ceff << " " << counter << " " << dLK << " " << dLK_sq << " " << Wr << " " << Wr_sq << " " << Tw << " " << Tw_sq << " " << L << " " << Wr_ful << " " << Wr_sq_ful << std::endl;
    }
    else {
        ofstr << chain->get_force() << " " << Ceff << " " << counter << " " << dLK << " " << dLK_sq << " " << Wr << " " << Wr_sq << " " << Tw << " " << Tw_sq << " " << L << std::endl;
//    ofstr << chain->get_force() << " " << Ceff << endl;
    }
    ofstr.close();

}

