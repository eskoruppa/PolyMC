#include "Dump_PersistenceLength.h"

Dump_PersLen::Dump_PersLen(Chain * ch, int N_dump, const std::string& filename, int max_dist,bool cal_tors, bool append)
: Dump(ch,N_dump,filename,append), m_max(max_dist), cal_torsional_perslen(cal_tors)
{
    if (max_dist < num_bp) m_max = max_dist;
    else                   m_max = num_bp-1;
//    if (max_dist < num_bp-1) m_max = max_dist;
//    else                   m_max = num_bp-2;
    tangent_dots       = arma::zeros(m_max);
    tangent_dots_count = arma::zeros(m_max);

    if (cal_torsional_perslen) {
        cosom3_sums       = arma::zeros(m_max);
        cosom3_sums_count = arma::zeros(m_max);
    }



}
Dump_PersLen::~Dump_PersLen() {}


//void Dump_PersLen::prod_dump() {
//    double val;
//    for (int n=0;n<m_max;n++) {
//        for (int i=0;i<(num_bp-m_max-1);i++) {
//            val = arma::dot(chain->get_triads()->slice(i).col(2),chain->get_triads()->slice(i+n+1).col(2));
//            if (val==val) {
//                tangent_dots(n)       += val;
//                tangent_dots_count(n) += 1;
//            }
//        }
//    }
//}

void Dump_PersLen::prod_dump() {
    double val;
    for (int n1=0;n1<num_bp-1;n1++) {
        int mlast = m_max;
        if (n1+m_max>=num_bp) {
            mlast = num_bp-1-n1;
        }
        for (int m=1;m<=mlast;m++) {
            val = arma::dot(chain->get_triads()->slice(n1).col(2),chain->get_triads()->slice(n1+m).col(2));
            tangent_dots(m-1)       += val;
            tangent_dots_count(m-1) += 1;
        }
    }
    if (cal_torsional_perslen) {
        double sum = 0;
        arma::colvec Theta;
        for (int n1=0;n1<num_bp-1;n1++) {
            sum = 0;
            int mlast = m_max;
            if (n1+m_max>=num_bp) {
                mlast = num_bp-1-n1;
            }
            for (int m=0;m<mlast;m++) {
                Theta = *BPS[n1+m]->get_Theta();
                sum += Theta(2);

                cosom3_sums(m)       += std::cos(sum);
                cosom3_sums_count(m) += 1;
            }
        }
    }
}

void Dump_PersLen::final_dump() {

	std::ofstream ofstr;
	if (app) ofstr.open(fn, std::ofstream::out | std::ofstream::app);
	else     ofstr.open(fn, std::ofstream::out | std::ofstream::trunc);

	arma::colvec lb = arma::zeros(m_max);
	for (int n=0;n<m_max;n++) {
		lb(n) = -(n+1.0)*chain->get_disc_len()/(std::log( tangent_dots(n)/tangent_dots_count(n) ));
	}

    ofstr << m_max;
    for (int n=0;n<m_max;n++) {
        ofstr << " " << lb(n);
//        ofstr << " " << tangent_dots(n)/tangent_dots_count(n);
    }
    ofstr << std::endl;

    if (cal_torsional_perslen) {
        arma::colvec lt = arma::zeros(m_max);
        for (int n=0;n<m_max;n++) {
            lt(n) = -(n+1.0)*chain->get_disc_len() / (std::log(cosom3_sums(n)/cosom3_sums_count(n)));
        }
        ofstr << m_max;
        for (int n=0;n<m_max;n++) {
            ofstr << " " << lt(n);
        }
        ofstr << std::endl;
    }
    ofstr.close();
}

