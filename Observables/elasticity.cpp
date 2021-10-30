#include "../Chain.h"

double Chain::rel_extension(arma::colvec& z_dir) {
    return extension(z_dir)/(disc_len*num_bps);
}
double Chain::extension(arma::colvec& z_dir) {
    return arma::dot(z_dir, bp_pos.col(num_bp-1)-bp_pos.col(0));
}


arma::colvec Chain::cal_static_persistence_length(int from, int to, int m_max) {
    int N = to-from;
    int M = 1000*N;
    if (m_max > M-1) {
        m_max = M-1;
    }

    arma::cube stat_triads = arma::zeros(3,3,M);
    stat_triads.slice(0)   = arma::eye<arma::mat>(3,3);

    for (int i=1;i<M;i++) {
        stat_triads.slice(i) =  stat_triads.slice(i-1) * *BPS[i%N]->get_R0();
    }

    arma::colvec lbstat = arma::zeros(m_max);

    for (int m=0;m<=m_max-1;m++) {
        double val = 0;
        double counter = 0;
        for (int i=0;i<M-m_max;i++) {
            val += arma::dot(stat_triads.slice(i).col(2),stat_triads.slice(i+m+1).col(2));
            counter++;
        }
//        lbstat(m) = val/counter;
        lbstat(m) = -(m+1)*disc_len/ std::log(val/counter);
    }
    return lbstat;
}
