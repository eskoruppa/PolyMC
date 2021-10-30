#include "Dump_SupercoilingPhase.h"

Dump_SupercoilingPhase::Dump_SupercoilingPhase(Chain * ch, int N_dump, const std::string& filename, bool append)
: Dump(ch,N_dump,filename,append)
{
    if (!app) {
        std::ofstream ofstr;
        ofstr.open(fn, std::ofstream::out | std::ofstream::trunc);
        ofstr.close();
    }
}
Dump_SupercoilingPhase::~Dump_SupercoilingPhase() {}

void Dump_SupercoilingPhase::prod_dump() {
    std::ofstream file;
    file.open(fn, std::ofstream::app);

    arma::mat writhemat;
    chain->langowski_writhe_elements(&writhemat,seg_size);

    int N = writhemat.n_cols;

    int first_row = 2;
    int last_row  = N-3;
    int first_col = 2;
    int last_col  = N-3;

    writhemat = writhemat.submat( first_row, first_col, last_row, last_col );
    writhemat = neg2zero(writhemat);
    arma::colvec wrdens = cal_writhe_dens(writhemat);

    double total = arma::accu(writhemat);
    double wrloc = wr_locality(wrdens, 0.5, total);

    file << total << " " << wrloc << "\n";
    file.close();
}

void Dump_SupercoilingPhase::final_dump() {
    prod_dump();
}


arma::colvec Dump_SupercoilingPhase::cal_writhe_dens(const arma::mat & writhemat) {
    int N = writhemat.n_cols;
    arma::colvec wd = arma::zeros<arma::colvec>(2*N-3);

    double wrsum;
    for (int i=0;i<N-1;i++) {
        wrsum = 0;
        for (int j=0;j<i+1;j++) {
            wrsum += writhemat(j,i+1-j);
        }
        wd(i) = wrsum;
    }
    for (int i=1;i<N-1;i++) {
        wrsum = 0;
        for (int j=0;j<N-i;j++) {
            wrsum += writhemat(i+j,N-1-j);
        }
        wd(N-2+i) = wrsum;
    }

    arma::colvec wdred = arma::zeros<arma::colvec>(N-2);
    for (int i=0;i<N-2;i++) {
        wdred(i) = wd(2*i) + wd(2*i+1);
    }
    return wdred;
}


double Dump_SupercoilingPhase::wr_locality(arma::vec wrdens, double perc, double total) {

    int N = wrdens.n_elem;
    double mean = arma::mean(wrdens);
    if (mean < 0) {
        wrdens = - wrdens;
        total  = -total;
    }
    for (int i=0;i<N;i++) {
        if (wrdens(i) < 0) {
            wrdens(i) = 0;
        }
    }

    double accum = 0;
    double frac  = 0;
    int argm;
    for (int i=0;i<N;i++) {
        argm         = argmax(wrdens);
        accum       += wrdens(argm);
        wrdens(argm) = 0;
        if (accum*1./total >= perc) {
            frac = (i+1.)/N;
            break;
        }
    }
    if (frac == 0) {
        frac = 1;
    }
    std::cout << "Writhe frac = " <<  frac << std::endl;
    return frac;
}

int Dump_SupercoilingPhase::argmax(const arma::colvec& vec) {
    int N = vec.n_elem;
    int    max_arg = 0;
    double max_val = vec(0);

    for (int i=1;i<N;i++) {
        if (vec(i) > max_val) {
            max_arg = i;
            max_val = vec(i);
        }
    }
    return max_arg;
}

arma::mat Dump_SupercoilingPhase::neg2zero(const arma::mat& WM) {
    int N = WM.n_cols;
    arma::mat zWM = arma::zeros<arma::mat>(N,N);

    for (int i=1;i<N;i++) {
        for (int j=1;j<N;j++) {
            if (WM(i,j) >= 0) {
                zWM(i,j) = WM(i,j);
            }
        }
    }
    return zWM;
}


