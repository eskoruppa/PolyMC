#include "ProximityPlecFinder.h"

ProximityPlecFinder::ProximityPlecFinder() :
autodetect_thresholds(true)
{
}

ProximityPlecFinder::ProximityPlecFinder(double threshold_dist,int min_seg_dist,double frac_closest) :
threshold_dist(threshold_dist),
min_seg_dist(min_seg_dist),
frac_closest(frac_closest),
autodetect_thresholds(false)
{
}

ProximityPlecFinder::~ProximityPlecFinder() {
}


double ProximityPlecFinder::get_xp(const arma::mat * pos) {
//    std::vector<int> plecs = find_plecs(pos);
    std::vector<int> plecs = find_plecs_Neuman(pos);
    int plec_segs = 0;

    std::cout << "Neuman" << std::endl;
    for (int i=0;i<plecs.size();i+=2) {
        std::cout << " " << plecs[i] << " " << plecs[i+1] << std::endl;
        plec_segs += (plecs[i+1] - plecs[i]);
    }
    double xp = 1.*plec_segs/pos->n_cols;
    return xp;
}


std::vector<int> ProximityPlecFinder::find_plecs(const arma::mat * pos) {
    cal_closest(pos);

    std::cout << "########################" << std::endl;

    std::vector<int> plecs;
    int num_bp = pos->n_cols;
    int id = 0;
    int entrance,exit;
    while (id < num_bp-min_seg_dist) {
//        std::cout << "id = " << id << std::endl;
        if (closest_dist[id] < threshold_dist && id < closest_id[id]) {
            entrance = id;
            exit     = closest_id[id];

            std::cout << entrance << " " << exit << std::endl;

            int num_within = 0;
            int num_count  = 0;
            for (int pid=entrance+1;pid<exit;pid++) {
                std::cout << "closest_id " << pid << " - " << closest_id[pid] << " - " << closest_dist[pid] << std::endl;
                if (closest_id[pid] < 0) {
                    continue;
                }
                num_count++;
                if ((entrance <= closest_id[pid]) && (closest_id[pid] <= exit)) {
                    num_within++;
                }
            }
            double frac_within = 1.*num_within/num_count;
            std::cout << " frac = " << frac_within << std::endl;
            std::cout << num_within << std::endl;
            std::cout << num_count << std::endl;
            std::cout << frac_closest << std::endl;

            if (frac_within >= frac_closest) {
                std::cout << entrance << " " << exit << std::endl;
                plecs.push_back(entrance);
                plecs.push_back(exit);
            }
            id = exit;
        }
        id++;
    }
    return plecs;
}


void ProximityPlecFinder::cal_closest(const arma::mat * pos) {
    unsigned num_bp = pos->n_cols;

    double dist;
    double min_dist;
    int    min_id;

    double initdist = arma::norm(pos->col(1)-pos->col(0))*num_bp*1000;
    closest_id   = -1*arma::ones<arma::ivec>(num_bp);
    closest_dist = initdist*arma::ones<arma::vec>(num_bp);

    for (int id1=0;id1<num_bp-min_seg_dist;id1++) {
        for (int id2=id1+min_seg_dist;id2<num_bp;id2++) {
            dist = arma::norm(pos->col(id1)-pos->col(id2));
            if (dist < closest_dist[id1]) {
                closest_dist[id1] = dist;
                closest_id[id1]   = id2;
            }
            if (dist < closest_dist[id2]) {
                closest_dist[id2] = dist;
                closest_id[id2]   = id1;
            }
        }
    }
}



std::vector<int> ProximityPlecFinder::find_plecs_Neuman(const arma::mat * pos) {
    std::vector<int> plecs;
    int num_bp = pos->n_cols;
    int id1 = 0;
    double closest_dist,dist;
    int    closest_id;
    while (id1 < num_bp-min_seg_dist) {
        closest_dist = arma::norm(pos->col(id1+min_seg_dist)-pos->col(id1));
        closest_id   = id1+min_seg_dist;
        for (int id2=id1+min_seg_dist+1;id2<num_bp;id2++) {
            dist = arma::norm(pos->col(id2)-pos->col(id1));
            if (dist < closest_dist) {
                closest_dist = dist;
                closest_id   = id2;
            }
        }
        if (closest_dist < threshold_dist) {
            plecs.push_back(id1);
            plecs.push_back(closest_id);
            id1 = closest_id;
        }
        id1++;
    }
    return plecs;
}












