#include "Dump_SolenoidCorrelationFunction.h"

Dump_SolCor::Dump_SolCor(Chain * ch, int N_dump, const std::string& filename, arma::colvec z_dir, int k_max, int n_max, bool append)
: Dump(ch,N_dump,filename,append), zdir(z_dir), nmax(n_max), kmax(k_max)
{
    if (kmax*4 >= num_bp) kmax = (int)(num_bp*1./4);
    cors       = arma::zeros(kmax,nmax);
    cors_count = arma::zeros(kmax,nmax);
    zdir = zdir/arma::norm(zdir);
}
Dump_SolCor::~Dump_SolCor() {}


void Dump_SolCor::prod_dump() {

    if (use_tangents) {
        int loc_nmax;
        double C;
        int k4;
        double zdot;
        double val;
        arma::colvec a,b,t;
        for (int k=1;k<=kmax;k++) {
            k4 = 4*k;
            for (int i=1;i<num_bp-1-k4;i++) {
                loc_nmax = smaller((num_bp-i-1)/k4,nmax);
                t    = chain->get_triads()->slice(i).col(2);
                zdot = arma::dot(t,zdir);
                if (abs(zdot-1)>1e-6) {

                    a = t - zdot*zdir;
                    a = a/arma::norm(a);
                    b = arma::cross(a,zdir);
                    for (int n=0;n<loc_nmax;n++) {
        //                cout << "#########" << endl;
        //                cout << i+(4*n+1)*k << " " << i+(4*n+3)*k << endl;
        //                cout << k-1 << " " << n << endl;
        //                cout << kmax << " " << nmax << endl;
//
//                        val = arma::dot(b,chain->get_triads()->slice(i+(4*n+1)*k).col(2)) - arma::dot(b,chain->get_triads()->slice(i+(4*n+3)*k).col(2));

//                        val = chain->get_triads()->slice(i)(0,2)*chain->get_triads()->slice(i+k)(1,2);
                        val = chain->get_triads()->slice(i)(0,2)*chain->get_triads()->slice(i+k)(1,2) - chain->get_triads()->slice(i)(1,2)*chain->get_triads()->slice(i+k)(0,2);

                        cors(k-1,n)  += val;
                        if (val!=val) {
                            std::cout << "Nan detected! " << std::endl;
                            std::cout << val << std::endl;
                            std::cout << b.t();
                            std::cout << chain->get_triads()->slice(i+(4*n+1)*k).col(2).t();
                            std::cout << t.t();
                            std::cout << a.t();
                        }

                        cors_count(k-1,n) += 1;
                    }
                }
            }
        }
    }
    else {
        arma::mat ktan;
        int       n_ktan;
        int       loc_nmax;
        double    val;

        int k4;
        for (int k=1;k<=kmax;k++) {
            k4 = 4*k;
            n_ktan = num_bp-k;
            ktan   = arma::zeros(3,n_ktan);
            for (int i=0;i<n_ktan;i++) {
                ktan.col(i) = chain->get_bp_pos()->col(i+k) - chain->get_bp_pos()->col(i);
            }
            for (int i=0;i<n_ktan;i++) {
                loc_nmax = smaller((n_ktan-i-1)/4,nmax);
                for (int n=0;n<loc_nmax;n++) {
                    val = arma::dot(arma::cross(ktan.col(i),ktan.col(i+1+n*4)),zdir);
//                    val = ktan(0,i)*ktan(1,i+1+n*4);
//                    val = ktan(0,i)*ktan(1,i+k);
                    cors(k-1,n)  += val;
                    cors_count(k-1,n) += 1;
                }
            }
        }
    }
}


void Dump_SolCor::final_dump() {
	std::ofstream ofstr;
	if (app) ofstr.open(fn, std::ofstream::out | std::ofstream::app);
	else     ofstr.open(fn, std::ofstream::out | std::ofstream::trunc);

    for (int k=0;k<kmax;k++) {
        ofstr << k+1;
        for (int n=0;n<nmax;n++) {
            ofstr << " " << cors(k,n)/cors_count(k,n);
        }
        ofstr << std::endl;
    }
    ofstr.close();
}

