#include "Dump_ClosureHist.h"

Dump_ClosureHist::Dump_ClosureHist(        
        Chain * ch, 
        int N_dump, 
        const std::string& filename, 
        int num_dist_bins, 
        int num_angle_bins,
        int num_twist_bins,
        double dist_lower, 
        double dist_upper, 
        double angle_lower, 
        double angle_upper,
        double twist_lower, 
        double twist_upper,
        bool bin_costheta)
: Dump(ch,N_dump,filename,false), 
active(false),
dist_active(false),
angle_active(false),
twist_active(false),
num_active(0),
num_dist_bins(num_dist_bins), 
num_angle_bins(num_angle_bins), 
num_twist_bins(num_twist_bins), 
dist_lower(dist_lower), 
dist_upper(dist_upper), 
angle_lower(angle_lower), 
angle_upper(angle_upper), 
twist_lower(twist_lower), 
twist_upper(twist_upper), 
bin_costheta(bin_costheta)
{
    if (num_dist_bins > 0) {
        dist_active = true;
        num_active++;
        active = true;
        dist_range = dist_upper - dist_lower;
        if (dist_range <= 0) {
            throw std::invalid_argument("dist_lower is larger or equal to dist_upper.");
        }
    }
    if (num_angle_bins > 0) {
        angle_active = true;
        num_active++;
        active = true;
        angle_range = angle_upper - angle_lower;
        if (angle_range <= 0) {
            throw std::invalid_argument("angle_lower is larger or equal to angle_upper.");
        }
    }
    if (num_twist_bins > 0) {
        twist_active = true;
        num_active++;
        active = true;
        twist_range = twist_upper - twist_lower;
        if (twist_range <= 0) {
            throw std::invalid_argument("twist_lower is larger or equal to twist_upper.");
        }
    }

    // combination ids
    // 100 -> 0
    // 010 -> 1
    // 001 -> 2
    // 110 -> 3
    // 101 -> 4
    // 011 -> 5
    // 111 -> 6
    if (dist_active && !angle_active && !twist_active) {
        combination_id = 0;
        num_bins_1 = num_dist_bins;
        hist = arma::zeros<arma::ivec>(num_dist_bins);
    }
    if (!dist_active && angle_active && !twist_active) {
        combination_id = 1;
        num_bins_1 = num_angle_bins;
        hist = arma::zeros<arma::ivec>(num_angle_bins);
    }
    if (!dist_active && !angle_active && twist_active) {
        combination_id = 2;
        num_bins_1 = num_twist_bins;
        hist = arma::zeros<arma::ivec>(num_twist_bins);
    }
    if (dist_active && angle_active && !twist_active) {
        combination_id = 3;
        num_bins_1 = num_dist_bins;
        num_bins_2 = num_angle_bins;
        mathist = arma::zeros<arma::imat>(num_dist_bins,num_angle_bins);
    }
    if (dist_active && !angle_active && twist_active) {
        combination_id = 4;
        num_bins_1 = num_dist_bins;
        num_bins_2 = num_twist_bins;
        mathist = arma::zeros<arma::imat>(num_dist_bins,num_twist_bins);
    }
    if (!dist_active && angle_active && twist_active) {
        combination_id = 5;
        num_bins_1 = num_angle_bins;
        num_bins_2 = num_twist_bins;
        mathist = arma::zeros<arma::imat>(num_angle_bins,num_twist_bins);
    }
    if (dist_active && angle_active && twist_active) {
        combination_id = 6;
        num_bins_1 = num_dist_bins;
        num_bins_2 = num_angle_bins;
        num_bins_3 = num_twist_bins;
        cubehist = arma::zeros<arma::icube>(num_dist_bins,num_angle_bins,num_twist_bins);
    }
    invalid_counter = 0;
}
Dump_ClosureHist::~Dump_ClosureHist() {}

void Dump_ClosureHist::prod_dump() {

    if (!active) {return;}

    int dist_id, angle_id, twist_id;

    if (dist_active) {
        double dist = arma::norm(chain->get_bp_pos()->col(num_bp-1)-chain->get_bp_pos()->col(0));
        dist_id = (int)(num_dist_bins*(dist-dist_lower)/dist_range);
        if (dist_id < 0 || dist_id >= num_dist_bins) {
            invalid_counter++;
            return;
        }
    }
    if (angle_active) {
        double costheta = arma::dot(chain->get_triads()->slice(0).col(2),chain->get_triads()->slice(num_bp-1).col(2));
        angle_id;
        if (bin_costheta) {
            angle_id = (int)(num_angle_bins*(costheta-angle_lower)/angle_range);
        }
        else{
            angle_id = (int)(num_angle_bins*(std::acos(costheta)-angle_lower)/angle_range);
        }
        if (angle_id < 0 || angle_id >= num_angle_bins) {
            invalid_counter++;
            return;
        }
    }
    if (twist_active) {
        arma::vec Omega = ExtractTheta(chain->get_triads()->slice(num_bp-1).t() * chain->get_triads()->slice(0) );
        double twist = Omega(2);
        twist_id = (int)(num_twist_bins*(twist-twist_lower)/twist_range);
        if (twist_id < 0 || twist_id >= num_twist_bins) {
            invalid_counter++;
            return;
        }
    }

    // Only distance
    if (combination_id==0) {
        hist(dist_id)++;
        return;
    }
    // Only angle
    if (combination_id==1) {
        hist(angle_id)++;
        return;
    }
    // Only twist
    if (combination_id==2) {
        hist(twist_id)++;
        return;
    }
    // dist and angle
    if (combination_id==3) {
        mathist(dist_id,angle_id)++;
        return;
    }
    // dist and twist
    if (combination_id==4) {
        mathist(dist_id,twist_id)++;
        return;
    }
    // angle and twist
    if (combination_id==5) {
        mathist(angle_id,twist_id)++;
        return;
    }
    // dist, angle, and twist
    if (combination_id==6) {
        cubehist(dist_id,angle_id,twist_id)++;
        return;
    }
}

void Dump_ClosureHist::final_dump() {

    if (!active) {
        return;
    }
    ///////////////////////////////////////////////
    // Write histogram file
    if (num_active == 1) {
        std::ofstream histfile;
        histfile.open(fn+".closure_hist", std::ofstream::out | std::ofstream::trunc);
        histfile << hist(0);
        for (unsigned i=1;i<num_bins_1;i++){
            histfile << " " << hist(i);
        }
        histfile << "\n";
    } 
    else if (num_active == 2)
    {
        std::ofstream histfile;
        histfile.open(fn+".closure_hist", std::ofstream::out | std::ofstream::trunc);
        for (unsigned i=0;i<num_bins_1;i++){
            histfile << mathist(i,0);
            for (unsigned j=1;j<num_bins_2;j++){
                histfile << " " << mathist(i,j);
            }
            histfile << "\n";
        }
    }
    else {
        std::ofstream histfile;
        histfile.open(fn+".closure_hist", std::ofstream::out | std::ofstream::trunc);
        for (unsigned i=0;i<num_bins_1;i++){
            for (unsigned j=0;j<num_bins_2;j++){
                histfile << cubehist(i,j,0);
                for (unsigned k=1;k<num_bins_3;k++){
                    histfile << " " << cubehist(i,j,k);
                }
                histfile << "\n";
            }
        }
    }

    ///////////////////////////////////////////////
    // Write shape file
    std::ofstream shapefile;
    shapefile.open(fn+".closure_hist_shape", std::ofstream::out | std::ofstream::trunc);
    if (dist_active) {
        shapefile << num_dist_bins << "\n";
    }
    if (angle_active) {
        shapefile << num_angle_bins << "\n";
    }
    if (twist_active) {
        shapefile << num_twist_bins << "\n";
    }

    ///////////////////////////////////////////////
    // Write active meta data
    std::ofstream activefile;
    activefile.open(fn+".closure_hist_active", std::ofstream::out | std::ofstream::trunc);
    if (dist_active) {
        activefile << "dist " << num_dist_bins << "\n";
    }
    if (angle_active) {
        activefile << "angle " << num_angle_bins << "\n";
    }
    if (twist_active) {
        activefile << "twist " << num_twist_bins << "\n";
    }

    ///////////////////////////////////////////////
    // Write bin edges file
    std::ofstream edgefile;
    edgefile.open(fn+".closure_hist_edges", std::ofstream::out | std::ofstream::trunc);
    double binsize;
    if (dist_active) {
        binsize = dist_range / num_dist_bins;
        edgefile << dist_lower;
        for (int i=1;i<=num_dist_bins;i++) {
            edgefile << " " << dist_lower + i*binsize;
        }
        edgefile << "\n";
    }
    if (angle_active) {
        binsize = angle_range / num_angle_bins;
        edgefile << angle_lower;
        for (int i=1;i<=num_angle_bins;i++) {
            edgefile << " " << angle_lower + i*binsize;
        }
        edgefile << "\n";
    }
    if (twist_active) {
        binsize = twist_range / num_twist_bins;
        edgefile << twist_lower;
        for (int i=1;i<=num_twist_bins;i++) {
            edgefile << " " << twist_lower + i*binsize;
        }
        edgefile << "\n";
    }
    edgefile.close();

    // ///////////////////////////////////////////////
    // // Write histogram file
    // std::ofstream histfile;
    // histfile.open(fn+".closure_hist", std::ofstream::out | std::ofstream::trunc);
    // if (!angle_active) { 
    //     histfile << hist(0);
    //     for (unsigned i=1;i<num_dist_bins;i++){
    //         histfile << " " << hist(i);
    //     }
    //     histfile << "\n";
    // }
    // else {
    //     for (unsigned i=0;i<num_dist_bins;i++){
    //         histfile << mathist(i,0);
    //         for (unsigned j=1;j<num_angle_bins;j++){
    //             histfile << " " << mathist(i,j);
    //         }
    //         histfile << "\n";
    //     }
    // }
    // histfile.close();

    // ///////////////////////////////////////////////
    // // Write bin edges file
    // std::ofstream edgefile;
    // edgefile.open(fn+".closure_histedges", std::ofstream::out | std::ofstream::trunc);
    // double binsize = dist_range / num_dist_bins;
    // edgefile << dist_lower;
    // for (int i=1;i<=num_dist_bins;i++) {
    //     edgefile << " " << dist_lower + i*binsize;
    // }
    // edgefile << "\n";
    // if (angle_active) {
    //     binsize = angle_range / num_angle_bins;
    //     edgefile << angle_lower;
    //     for (int i=1;i<=num_angle_bins;i++) {
    //         edgefile << " " << angle_lower + i*binsize;
    //     }
    //     edgefile << "\n";
    // } 
    // edgefile.close();
}


bool Dump_ClosureHist::is_active() {
    return active;
}

