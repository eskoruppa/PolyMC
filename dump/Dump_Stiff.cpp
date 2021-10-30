#include "Dump_Stiff.h"

Dump_Stiff::Dump_Stiff(Chain * ch, int N_dump, const std::string& filename, bool all_bps)
: Dump(ch,N_dump,filename,false), counter(0), all(all_bps)
{
    cov  = arma::zeros(3,3);
    if (all_bps) covs = arma::zeros(3,3,num_bps);
}
Dump_Stiff::~Dump_Stiff() {}


void Dump_Stiff::prod_dump() {
    counter++;
    arma::mat C;
    for (int i=0;i<num_bps;i++) {
        C = *BPS[i]->get_Theta()*BPS[i]->get_Theta()->t();
        cov += C;
        if (all) {
            covs.slice(i) += C;
        }
    }
}

void Dump_Stiff::final_dump() {

    if (counter > 0) {
        std::ofstream ofstr;
        ofstr.open(fn, std::ofstream::out | std::ofstream::trunc);
        arma::mat Mat,Cov;

        ofstr << "MEAN" << std::endl;
        Cov = cov/(counter*num_bps);
        Mat = chain->get_disc_len()*arma::inv_sympd(Cov);

//        ofstr << "cov " << Cov(0,0) << " " << Cov(0,1) << " " << Cov(0,2) << " "
//                        << Cov(1,0) << " " << Cov(1,1) << " " << Cov(1,2) << " "
//                        << Cov(2,0) << " " << Cov(2,1) << " " << Cov(2,2) << endl;
//        ofstr << "mat " << Mat(0,0) << " " << Mat(0,1) << " " << Mat(0,2) << " "
//                        << Mat(1,0) << " " << Mat(1,1) << " " << Mat(1,2) << " "
//                        << Mat(2,0) << " " << Mat(2,1) << " " << Mat(2,2) << endl;
        ofstr << "cov" << std::endl;
        ofstr << Cov;
        ofstr << "mat" << std::endl;
        ofstr << Mat;

        std::cout << Mat;

        if (all) {
            ofstr << "\nALL" << std::endl;
            arma::mat Cov;
            for (int i=0;i<num_bps;i++) {
                Cov = covs.slice(i)/counter;
                Mat = chain->get_disc_len()*arma::inv_sympd(Cov);

                ofstr << "#" << i << std::endl;
//                ofstr << "cov " << Cov(0,0) << " " << Cov(0,1) << " " << Cov(0,2) << " "
//                                << Cov(1,0) << " " << Cov(1,1) << " " << Cov(1,2) << " "
//                                << Cov(2,0) << " " << Cov(2,1) << " " << Cov(2,2) << endl;
//                ofstr << "mat " << Mat(0,0) << " " << Mat(0,1) << " " << Mat(0,2) << " "
//                                << Mat(1,0) << " " << Mat(1,1) << " " << Mat(1,2) << " "
//                                << Mat(2,0) << " " << Mat(2,1) << " " << Mat(2,2) << endl;
                ofstr << "cov" << std::endl;
                ofstr << Cov;
                ofstr << "mat" << std::endl;
                ofstr << Mat;
            }
        }
        ofstr.close();
    }
}

