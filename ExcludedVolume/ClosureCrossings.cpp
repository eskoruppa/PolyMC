#include "ClosureCrossings.h"

double front_closure(   const arma::colvec p1,
                        const arma::colvec p2,
                        const arma::colvec r0,
                        const arma::colvec y,
                        double tol) {

    arma::colvec v = p2-p1;

    /*
        Check if p1 and p2 are the same
    */
    double vsq = arma::dot(v,v);
    if (vsq < tol) {
        return distance_point_semi_infinite_lineseg(p1, r0, y);
    }

    /*
        Check if vectors are parallel
    */
    if ( vectors_parallel(v,y)) {
        return distance_parallel_linesegs_finite_semiinf(p1, p2,r0, y);
    }

    /*
        calculate lambda0 and mu0
    */

    double vy      = arma::dot(v,y);
    arma::colvec a = v - vy*y;
    arma::colvec b = p1-r0;
    double ay      = arma::dot(a,y);

    double lam0 = ( arma::dot(a,b)-arma::dot(b,y)*ay ) / ( vy*ay-arma::dot(a,v) );
    double mu0 = arma::dot(b + lam0*v,y);

    /*
        Evaluate different cases
    */
    if (mu0 < 0) {
        return distance_point_lineseg(r0,p1,p2);
    }

    if (lam0 < 0) {
        return distance_point_semi_infinite_lineseg(p1,r0,y);
    }

    if (lam0 > 1) {
        return distance_point_semi_infinite_lineseg(p2,r0,y);
    }

    arma::colvec plam0 = p1 + lam0*v;
    arma::colvec lmu0  = r0 + mu0*y;
    return arma::norm(plam0-lmu0);
}

double back_closure(const arma::colvec p1,
                    const arma::colvec p2,
                    const arma::colvec r_l1,
                    const arma::colvec r_l2,
                    const arma::colvec y,
                    double tol) {

    arma::colvec p2p = p2 + r_l1 - r_l2;
    return front_closure( p1, p2p, r_l1, y, tol);
}
