#include "LinkingNumber.h"

double fullerwrithe_vertex_angle(const arma::colvec& a,const arma::colvec& b,const arma::colvec& c) {
    arma::colvec cross1 = arma::cross(a,b);
    arma::colvec cross2 = arma::cross(a,c);
    double dot = arma::dot(normalize(cross1),normalize(cross2));
    double A = std::acos(dot);
    if (A!=A) {
        A = 0;
    }
    return A;
}

double fullerwrithe_triangle_area(const arma::colvec& ez,const arma::colvec& t1,const arma::colvec& t2) {
    double A = fullerwrithe_vertex_angle(ez,t1,t2);
    double B = fullerwrithe_vertex_angle(t1,ez,t2);
    double C = fullerwrithe_vertex_angle(t2,t1,ez);
    double abs_Omega = std::abs(A+B+C-M_PI);
    double sgn_Omega = sgn(arma::dot(arma::cross(ez,t1),t2));
    double Omega = abs_Omega*sgn_Omega;
    return Omega;
}

double fullerwrithe_writhesum(const arma::mat& tans,const arma::colvec& z_dir) {
    int N    = tans.n_cols;
    arma::colvec ez = normalize(z_dir);

    double wrsum = 0;
    for (unsigned i=0;i<N-1;i++){
        wrsum += fullerwrithe_triangle_area(ez,tans.col(i),tans.col(i+1));
    }
    double FWr = wrsum / (2*M_PI);
    return FWr;
}

double cal_closure_writhe(const arma::mat& pos,const arma::colvec& z_dir) {
    // contributions to the writhe stemming from the closure
    arma::colvec ez = normalize(z_dir);
    arma::colvec r0 = pos.col(0);

    arma::mat U = pos;
    for (unsigned i=0;i<pos.n_cols;i++){
        U.col(i) = normalize(U.col(i)-r0);
    }
    double closure_writhe = 0;
    closure_writhe += fullerwrithe_writhesum( U,ez);

    r0 = pos.col(pos.n_cols-1);
    for (unsigned i=0;i<pos.n_cols;i++){
        U.col(i) = -normalize(U.col(i)-r0);
    }
    closure_writhe += fullerwrithe_writhesum(-U,ez);
    return closure_writhe;
}

double cal_fullerwrithe(const arma::mat& pos,const arma::colvec& z_dir) {
////////////////////////////////////////////////////////////////////////////////////
// Calculation of fuller writhe using the method outlined in Chou 2014. Includes the
// contributions from the chain closure and calculates the fuller single integral in
// terms of the angle swept out by the tangents on the unit sphere.
////////////////////////////////////////////////////////////////////////////////////

    int Npos = pos.n_cols;
    int N    = Npos-1;
    arma::colvec ez = normalize(z_dir);

    // writhe of the center
    arma::mat tans = arma::zeros(3,N);
    arma::colvec tan;
    for (unsigned i=0;i<N;i++){
        tan = pos.col(i+1) - pos.col(i);
        tans.col(i) = normalize(tan);
    }
    double wrsum_hh       = fullerwrithe_writhesum(tans,ez);
    double closure_writhe = cal_closure_writhe(pos,ez);

//    double FWr = wrsum_hh + closure_writhe;
    double FWr = wrsum_hh;
    return FWr;
}



///////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////


double cal_langowskiwrithe(const arma::mat& pos,bool closed) {
    int N = pos.n_cols;
    double Omega_sum = 0;
    double Omega;
    for (unsigned i=0;i<N-3;i++) {
        for (unsigned j=i+2;j<N-1;j++) {
//            std::cout << i << " " << j << std::endl;
//            std::cout << pos.col(j).t() << std::endl;
//            std::cout << pos.col(j+1).t() << std::endl;
//            std::cout << pos.col(i).t() << std::endl;
//            std::cout << pos.col(i+1).t() << std::endl;

            Omega = langowskiwrithe_quadrangle_area(pos.col(i),pos.col(i+1),pos.col(j),pos.col(j+1));
//
//            if ( j == i + 1 ) {
//                std::cout << Omega << std::endl;
//            }

            Omega_sum += Omega;
        }
    }
    return Omega_sum/(2*M_PI);
}


double langowskiwrithe_quadrangle_area (const arma::colvec& r1,const arma::colvec& r2,const arma::colvec& r3,const arma::colvec& r4) {
    arma::colvec r12,r13,r14,r23,r24,r34;
    arma::colvec n1,n2,n3,n4;

    r12 = r2-r1;
    r13 = r3-r1;
    r14 = r4-r1;
    r23 = r3-r2;
    r24 = r4-r2;
    r34 = r4-r3;

    n1 = langowskiwrithe_nvec(r13,r14);
    n2 = langowskiwrithe_nvec(r14,r24);
    n3 = langowskiwrithe_nvec(r24,r23);
    n4 = langowskiwrithe_nvec(r23,r13);

    double alpha = langowskiwrithe_angle_contribution(n1,n2);
    double beta  = langowskiwrithe_angle_contribution(n2,n3);
    double gamma = langowskiwrithe_angle_contribution(n3,n4);
    double delta = langowskiwrithe_angle_contribution(n4,n1);

    double abs_Omega = std::abs(alpha+beta+gamma+delta);
    double sgn_Omega = sgn(arma::dot(arma::cross(r34,r12),r13));
    return abs_Omega*sgn_Omega;
}



arma::colvec langowskiwrithe_nvec (const arma::colvec& a,const arma::colvec& b) {
    arma::colvec n    = arma::cross(a,b);
    double       norm = arma::norm(n);
    if (std::abs(norm) < LINKINKNUMBER_NORMZERO_THRESHOLD) {
        return arma::zeros(3);
    }
    return n/norm;
}

double langowskiwrithe_angle_contribution (const arma::colvec& a,const arma::colvec& b) {
    double dot = arma::dot(a,b);
    if (dot > 1) {
        dot = 1;
    }
    else if (dot < -1) {
        dot = -1;
    }
    double angle = std::asin(dot);
    if (angle!=angle) {
        angle = 0;
    }
    return angle;
}


///////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////


double cal_langowskiwrithe_1b(const arma::mat& pos,bool closed) {
    int N = pos.n_cols;
    double Omega_sum = 0;
    double Omega;

//    std::cout << "--------------------" << std::endl;
//    std::cout << pos.col(N-2).t() << std::endl;
//    std::cout << pos.col(N-1).t() << std::endl;
//
//    std::cout << pos.col(0).t() << std::endl;
//    std::cout << pos.col(1).t() << std::endl;

    for (unsigned i=0;i<N-3;i++) {
        for (unsigned j=i+2;j<N-1;j++) {
//            std::cout << i << " " << j << std::endl;
//            std::cout << pos.col(j).t() << std::endl;
//            std::cout << pos.col(j+1).t() << std::endl;
//            std::cout << pos.col(i).t() << std::endl;
//            std::cout << pos.col(i+1).t() << std::endl;

            Omega = langowskiwrithe1b_Omega(pos.col(i),pos.col(i+1),pos.col(j),pos.col(j+1));
//
//            if ( j == i + 1 ) {
//                std::cout << Omega << std::endl;
//            }

            Omega_sum += Omega;
        }
    }
    return Omega_sum/(2*M_PI);
}

double langowskiwrithe1b_Omega(const arma::colvec& r1,const arma::colvec& r1p,const arma::colvec& r2,const arma::colvec& r2p) {
    arma::colvec vs1,vs2;
    arma::colvec e1,e2,r12;
    double s1,s2;
    double cosbeta,sinsqbeta;
    double a0,a1,a2;
    double Omega;

    vs1 = r1p-r1;
    vs2 = r2p-r2;
    s1 = arma::norm(vs1);
    s2 = arma::norm(vs2);

    e1 = vs1/s1;
    e2 = vs2/s2;

    r12 = r2-r1;

    cosbeta   = arma::dot(e1,e2);
    sinsqbeta = 1 - cosbeta*cosbeta;

    a1 = arma::dot(r12,cosbeta*e2-e1)/sinsqbeta;
    a2 = arma::dot(r12,e2-e1*cosbeta)/sinsqbeta;
    a0 = arma::dot(r12,arma::cross(e1,e2))/sinsqbeta;

    Omega  = langowskiwrithe1b_F(a1+s1,a2+s2,a0,cosbeta,sinsqbeta);
    Omega -= langowskiwrithe1b_F(a1+s1,a2,a0,cosbeta,sinsqbeta);
    Omega -= langowskiwrithe1b_F(a1,a2+s2,a0,cosbeta,sinsqbeta);
    Omega += langowskiwrithe1b_F(a1,a2,a0,cosbeta,sinsqbeta);

    return Omega;
}



double langowskiwrithe1b_F(double t1,double t2,double a0,double cosbeta,double sinsqbeta) {
    double nom   = t1*t2+a0*a0*cosbeta;
    double denom = a0*std::sqrt(t1*t1+t2*t2-2*t1*t2*cosbeta+a0*a0*sinsqbeta);
    double F     = -std::atan(nom/denom);
    if (F!=F) {
        return 0;
    }
    return F;
}

