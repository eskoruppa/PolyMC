#include "Dump_TangentCorrelation.h"

Dump_TanCor::Dump_TanCor(Chain * ch, int N_dump, const std::string& filename, int max_dist, bool append)
: Dump(ch,N_dump,filename,append), m_max(max_dist)
{
    if (max_dist < num_bp) m_max = max_dist;
    else                   m_max = num_bp-1;
    tangent_dots       = arma::zeros(m_max);
    tangent_dots_count = arma::zeros(m_max);
}
Dump_TanCor::~Dump_TanCor() {}


void Dump_TanCor::prod_dump() {
    for (int n=0;n<m_max;n++) {
        for (int i=0;i<(num_bp-m_max-1);i++) {
            tangent_dots(n)       += arma::dot(chain->get_triads()->slice(i).col(2),chain->get_triads()->slice(i+n+1).col(2));
            tangent_dots_count(n) += 1;
        }
    }
}


void Dump_TanCor::final_dump() {

    arma::colvec tancor = arma::zeros(m_max);
    for (int n=0;n<m_max;n++) {
        tancor(n) = tangent_dots(n)/tangent_dots_count(n);
    }

	std::ofstream ofstr;
	if (app) ofstr.open(fn, std::ofstream::out | std::ofstream::app);
	else     ofstr.open(fn, std::ofstream::out | std::ofstream::trunc);

    ofstr << m_max;
    for (int n=0;n<m_max;n++) {
        ofstr << " " << tancor(n);
    }
    ofstr << std::endl;
    ofstr.close();
}

