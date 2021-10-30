#include "Dump_Thetas.h"

Dump_Thetas::Dump_Thetas(Chain * ch, int N_dump, const std::string& filename,bool append)
: Dump(ch,N_dump,filename,append), counter(0)
{
    if (!app) {
        std::ofstream ofstr;
        ofstr.open(fn, std::ofstream::out | std::ofstream::trunc);
        ofstr.close();
    }
}
Dump_Thetas::~Dump_Thetas() {}

void Dump_Thetas::prod_dump() {
    counter++;

    arma::colvec Theta;
    std::ofstream ofstr;
    ofstr.open(fn, std::ofstream::out | std::ofstream::app);

//    ofstr << "#" << counter << " number_bps " << num_bps << "\n";
    for (int i=0;i<num_bps;i++) {
//        Theta = *BPS[i]->get_Theta();
        Theta = BPS[i]->get_FullTheta();
        ofstr << i << " " << Theta(0) << " " << Theta(1) << " " << Theta(2) << " " << "\n";
    }
    ofstr.close();
}


Dump_avgThetas::Dump_avgThetas(Chain * ch, int N_dump, const std::string& filename,bool append)
: Dump(ch,N_dump,filename,append), counter(0)
{
    Thetas      = arma::zeros(3,num_bps);
    Thetas_sq   = arma::zeros(3,num_bps);
    kappa       = arma::zeros(num_bps);
}
Dump_avgThetas::~Dump_avgThetas() {}

void Dump_avgThetas::prod_dump() {
    counter++;
    arma::colvec Theta;
    for (int i=0;i<num_bps;i++) {
        Theta = *BPS[i]->get_Theta();
        Thetas.col(i)    += Theta;
        Thetas_sq.col(i) += Theta%Theta;
        kappa(i)         += std::sqrt(Theta(0)*Theta(0)+Theta(1)*Theta(1));
    }
}


void Dump_avgThetas::final_dump() {

    std::ofstream ofstr;
    if (app) ofstr.open(fn, std::ofstream::out | std::ofstream::app);
    else     ofstr.open(fn, std::ofstream::out | std::ofstream::trunc);

    arma::colvec Theta;
    arma::colvec Thetasq;
//    ofstr << "N " << counter << " number_bp " << num_bps << "\n";
    for (int i=0;i<num_bps;i++) {
        Thetas.col(i)    = Thetas.col(i)/counter;
        kappa(i)         = kappa(i)/counter;
        Thetas_sq.col(i) = Thetas_sq.col(i)/counter;
        Theta   = Thetas.col(i);
        Thetasq = Thetas_sq.col(i);
        ofstr << i << " " << Theta(0) << " " << Theta(1) << " " << Theta(2) << " " << Thetasq(0) << " " << Thetasq(1) << " " << Thetasq(2) << " " << kappa(i) <<  "\n";
    }
    ofstr.close();
}

