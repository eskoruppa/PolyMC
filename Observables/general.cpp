#include "../Chain.h"


arma::mat Chain::distmap(double density) {
    int num_points = (int)(num_bp*density);
    int point_dist = num_bp/num_points;
    double dist;
    arma::mat distmat(num_points,num_points);
    for (int i=0;i<num_points-1;i++) {
        distmat(i,i) = 0;
        for (int j=i+1;j<num_points;j++) {
            dist = arma::norm(bp_pos.col(i*point_dist)-bp_pos.col(j*point_dist));
            distmat(i,j) = dist;
            distmat(j,i) = dist;
        }
    }
    return distmat;
}


