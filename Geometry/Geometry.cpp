#include "Geometry.h"

double distance_point_lineseg(const arma::colvec p, const arma::colvec q1, const arma::colvec q2) {
    arma::colvec v = q2 - q1;
    arma::colvec w = p  - q1;
    double lambda = arma::dot(w,v) / arma::dot(v,v);
    if (lambda < 0) {
        return arma::norm(q1-p);
    }
    if (lambda > 1) {
        return arma::norm(q2-p);
    }
    return arma::norm(q1+lambda*v-p);
}

double distance_point_semi_infinite_lineseg(const arma::colvec p, const arma::colvec q1, const arma::colvec v) {
    /*
        Distance between point p and semi-infinite linesegment bouned by the point q1 and extending in the direction v
    */
    arma::colvec w = p  - q1;
    double lambda = arma::dot(w,v) / arma::dot(v,v);
    if (lambda < 0) {
        return arma::norm(q1-p);
    }
    return arma::norm(q1+lambda*v-p);
}

bool vectors_parallel(const arma::colvec v1, const arma::colvec v2, double tol) {
    /*
        Check if vectors are parallel
    */
    double x1,x2,x3;
    if (std::abs(v2(0)) > 1e-10) {
        x1 = v1(0)/v2(0);
        x2 = v2(1)*x1;
        x3 = v2(2)*x1;
    }
    else if (std::abs(v2(1)) > 1e-10) {
        x1 = v1(1)/v2(1);
        x2 = v2(0)*x1;
        x3 = v2(2)*x1;
    }
    else if (std::abs(v2(2)) > 1e-10) {
        x1 = v1(2)/v2(2);
        x2 = v2(0)*x1;
        x3 = v2(1)*x1;
    }
    else return false;


    if (std::abs(x2-v1(1)) > tol) {
        return false;
    }
    if (std::abs(x2-v1(1)) > tol) {
        return false;
    }
    return true;
}

double distance_parallel_linesegs_finite_semiinf(const arma::colvec p1, const arma::colvec p2, const arma::colvec q, const arma::colvec v) {
    arma::colvec w1 = p1-q;
    arma::colvec w2 = p2-q;

    if (arma::dot(w1,v) < 0 && arma::dot(w2,v) < 0) {
        double dist1,dist2;
        dist1 = arma::norm(w1);
        dist2 = arma::norm(w2);
        if (dist1 < dist2) return dist1;
        return dist2;
    }
    arma::colvec uv = v / arma::norm(v);
    arma::colvec vd = w1 - arma::dot(w1,uv)*uv;
    return arma::norm(vd);
}
