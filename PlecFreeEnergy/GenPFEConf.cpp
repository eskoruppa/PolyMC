#include "GenPFEConf.h"


int gen_PFE_conf(arma::mat * pos, arma::cube * triads, int num_segs, double  disc_len, double termini_dist, double EV_rad, double loop_frac) {
    /*
        returns the index of the monomer at the boundary of the plectomene. Monomers of lower index constitute the bridge intended to
        preserve the linking number of the chain. These monomers should be excluded from any MC moves.
    */

    int N_con = 0;
    if (EV_rad > 0) {
        double Chi   = termini_dist+(1-1./std::sqrt(2))*disc_len;
        double discr = EV_rad*EV_rad-0.5*disc_len*disc_len;
        if (discr < 0) {
            throw std::logic_error("EV_rad smaller than disc_len");
        }
        N_con = int(std::floor((Chi-std::sqrt(discr))/disc_len));
    }

    double MAX_PHI = 0.75*M_PI;

    bool even = (num_segs%2==0);
    double  L = (num_segs-1)*disc_len;
    double half_tdist = 0.5*termini_dist;

    double threshold_gamma = (3*M_PI+4)*termini_dist/(4*L);
    double gamma;
    if (loop_frac < threshold_gamma) gamma = threshold_gamma;
    else                             gamma = loop_frac;

    double R = (2*gamma*L)/(3*M_PI+4);
    int left,right;

    arma::colvec point;
    arma::mat pts = arma::zeros(3,num_segs+N_con);

    // generate circle part
    double phi;
    double dphi = 2*std::abs(std::asin(disc_len/(2*R)));
    if (even) {
        phi  = 0.5*dphi;
        right = int(num_segs/2);
        left  = right-1;
    }
    else {
        phi = dphi;
        int mid = int((num_segs-1)/2);
        point = {0,R,0};
        pts.col(mid+N_con) = point;
        left  = mid-1;
        right = mid+1;
    }

    double x,y;
    while (phi < MAX_PHI) {
        x = R*std::sin(phi);
        y = R*std::cos(phi);
        if ((y < 0) && (x < half_tdist)) break;

        point = {-x,y,0};
        pts.col(left+N_con) = point;
        point = {x,y,0};
        pts.col(right+N_con) = point;

        left--;
        right++;
        phi += dphi;
    }

    // generate straight connectors
    double dx,dy,trsh;
    dx   = disc_len/std::sqrt(2);
    trsh = half_tdist+dx;
    while (pts(0,right-1+N_con) > trsh) {
        pts(0,left+N_con)  = pts(0,left+1+N_con) + dx;
        pts(1,left+N_con)  = pts(1,left+1+N_con) - dx;
        pts(0,right+N_con) = pts(0,right-1+N_con) - dx;
        pts(1,right+N_con) = pts(1,right-1+N_con) - dx;

        left-- ;
        right++;
    }

    dx = pts(0,right-1+N_con)-half_tdist;
    dy = std::sqrt(disc_len*disc_len-dx*dx);

    pts(0,left+N_con)  = -half_tdist;
    pts(1,left+N_con)  = pts(1,left+1+N_con) - dy;
    pts(0,right+N_con) =  half_tdist;
    pts(1,right+N_con) = pts(1,right-1+N_con) - dy;

    left--;
    right++;

    // generate straight parallel termini
    while (left >= 0) {
        pts(0,left+N_con)  = -half_tdist;
        pts(1,left+N_con)  = pts(1,left+1+N_con) - disc_len;
        pts(0,right+N_con) =  half_tdist;
        pts(1,right+N_con) = pts(1,right-1+N_con) - disc_len;

        left--;
        right++;
    }

    // generate terminus connectors
    if (N_con > 0) {
        double refx = pts(0,N_con) + disc_len/std::sqrt(2) - disc_len;
        double refy = pts(1,N_con) - disc_len/std::sqrt(2);
        for (int i=0;i<N_con;i++) {
            pts(0,i) = refx + (N_con-i)*disc_len;
            pts(1,i) = refy;
        }
    }

    // testing
    std::cout << "######" << std::endl;
    arma::colvec diff;
    for (int i=0;i<pts.n_cols-1;i++) {
        diff = pts.col(i+1)-pts.col(i);
        double n = arma::norm(diff);
        if (std::abs(n-disc_len)>1e-12) {
            std::cout << i << " " << n << std::endl;
        }
    }

    /*
        Calcuate Triads
    */
    arma::colvec e1,e2,e3;
    e2 = {0,0,1};

    arma::cube trds = arma::zeros(3,3,pts.n_cols);
    for (int i=0;i<pts.n_cols-1;i++) {
        e3 = pts.col(i+1)-pts.col(i);
        e3 = e3/arma::norm(e3);
        e1 = arma::cross(e2,e3);

        trds.slice(i).col(0) = e1;
        trds.slice(i).col(1) = e2;
        trds.slice(i).col(2) = e3;
    }
    e3 = {0,-1,0};
    e1 = arma::cross(e2,e3);
    trds.slice(pts.n_cols-1).col(0) = e1;
    trds.slice(pts.n_cols-1).col(1) = e2;
    trds.slice(pts.n_cols-1).col(2) = e3;

    *pos    = pts;
    *triads = trds;
    return N_con;
}
