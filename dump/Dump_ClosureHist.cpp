#include "Dump_ClosureHist.h"

Dump_ClosureHist::Dump_ClosureHist(        
        Chain * ch, 
        int N_dump, 
        const std::string& filename, 
        int num_dist_bins, 
        int num_angle_bins, 
        double dist_lower, 
        double dist_upper, 
        double angle_lower, 
        double angle_upper,
        bool bin_costheta)
: Dump(ch,N_dump,filename,false), angle_active(false), num_dist_bins(num_dist_bins), num_angle_bins(num_angle_bins), dist_lower(dist_lower), dist_upper(dist_upper), angle_lower(angle_lower), angle_upper(angle_upper), bin_costheta(bin_costheta)
{

    if (num_dist_bins <= 0) {
        active = false;
        return;
    }
    
    dist_range = dist_upper - dist_lower;
    angle_range = angle_upper - angle_lower;

    if (num_angle_bins > 0 && angle_range > 0) {
        angle_active = true;
        mathist = arma::zeros<arma::imat>(num_dist_bins,num_angle_bins);
    }
    else {
        if (dist_range < 0 || num_dist_bins <= 0) {
            active = false;
        }
        hist = arma::zeros<arma::ivec>(num_dist_bins);
    }
    invalid_counter = 0;
}
Dump_ClosureHist::~Dump_ClosureHist() {}

void Dump_ClosureHist::prod_dump() {

    if (!active) {return;}
    double dist = arma::norm(chain->get_bp_pos()->col(num_bp-1)-chain->get_bp_pos()->col(0));
    int dist_id, angle_id;
    dist_id = (int)(num_dist_bins*(dist-dist_lower)/dist_range);
    if (dist_id < 0 || dist_id >= num_dist_bins) {
        return;
    }
    if (!angle_active) {
        hist(dist_id)++;
        return;
    }
    double costheta = arma::dot(chain->get_triads()->slice(0).col(2),chain->get_triads()->slice(num_bp-1).col(2));
    if (bin_costheta) {
        angle_id = (int)(num_angle_bins*(costheta-angle_lower)/angle_range);
    }
    else{
        angle_id = (int)(num_angle_bins*(std::acos(costheta)-angle_lower)/angle_range);
    }
    if (angle_id < 0 || angle_id >= num_angle_bins) {
        return;
    }
    mathist(dist_id,angle_id)++;
}

void Dump_ClosureHist::final_dump() {

    if (!active) {
        return;
    }
    ///////////////////////////////////////////////
    // Write histogram file
    std::ofstream histfile;
    histfile.open(fn+".closure_hist", std::ofstream::out | std::ofstream::trunc);
    if (!angle_active) { 
        histfile << hist(0);
        for (unsigned i=1;i<num_dist_bins;i++){
            histfile << " " << hist(i);
        }
        histfile << "\n";
    }
    else {
        for (unsigned i=0;i<num_dist_bins;i++){
            histfile << mathist(i,0);
            for (unsigned j=0;j<num_angle_bins;j++){
                histfile << " " << mathist(i,j);
            }
            histfile << "\n";
        }
    }
    histfile.close();

    ///////////////////////////////////////////////
    // Write bin edges file
    std::ofstream edgefile;
    edgefile.open(fn+".closure_histedges", std::ofstream::out | std::ofstream::trunc);
    double binsize = dist_range / num_dist_bins;
    edgefile << dist_lower;
    for (int i=1;i<=num_dist_bins;i++) {
        edgefile << " " << dist_lower + i*binsize;
    }
    edgefile << "\n";
    if (angle_active) {
        binsize = angle_range / num_angle_bins;
        edgefile << angle_lower;
        for (int i=1;i<=num_angle_bins;i++) {
            edgefile << " " << angle_lower + i*binsize;
        }
        edgefile << "\n";
    } 
    edgefile.close();
}


bool Dump_ClosureHist::is_active() {
    return active;
}

