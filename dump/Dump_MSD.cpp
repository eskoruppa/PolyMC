#include "Dump_MSD.h"

Dump_MSD::Dump_MSD(Chain * ch, int N_dump, const std::string& filename, int max_dist, bool append)
: Dump(ch,N_dump,filename,append), m_max(max_dist)
{
    if (max_dist < num_bp) m_max = max_dist;
    else                   m_max = num_bp-1;
    sqdis       = arma::zeros(m_max);
    sqdis_count = arma::zeros(m_max);
}
Dump_MSD::~Dump_MSD() {}


void Dump_MSD::prod_dump() {
    double dist;
    for (int n=0;n<m_max;n++) {
        for (int i=0;i<(num_bp-m_max-1);i++) {
            dist = arma::norm( chain->get_bp_pos()->col(i+n+1) - chain->get_bp_pos()->col(i) );
            sqdis(n)       += dist*dist;
            sqdis_count(n) += 1;
        }
    }
}


void Dump_MSD::final_dump() {

	arma::colvec msd = arma::zeros(m_max);
	for (int n=0;n<m_max;n++) {
		msd(n) = sqdis(n)/sqdis_count(n);
	}

	std::ofstream ofstr;
	if (app) ofstr.open(fn, std::ofstream::out | std::ofstream::app);
	else     ofstr.open(fn, std::ofstream::out | std::ofstream::trunc);

    ofstr << m_max;
    for (int n=0;n<m_max;n++) {
        ofstr << " " << msd(n);
    }
    ofstr << std::endl;
    ofstr.close();
}

