#include "Dump_PlasmidEndpoints.h"

Dump_PlasmidEndpoints::Dump_PlasmidEndpoints(Chain * ch, int N_dump,
                                                         const std::string& filename,
                                                         double density,
                                                         int sample_step_dist,
                                                         int sample_points,
                                                         bool dump_writhe,
                                                         bool dump_xyz,
                                                         bool append)

: Dump(ch,N_dump,filename,append), density(density),sample_step_dist(sample_step_dist),sample_points(sample_points),dump_writhe(dump_writhe),dump_xyz(dump_xyz)
{
    if (!app) {

        std::ofstream ofstr;
        ofstr.open(fn, std::ofstream::out | std::ofstream::trunc);
        ofstr.close();

        ofstr.open(fn+".writhe", std::ofstream::out | std::ofstream::trunc);
        ofstr.close();

        ofstr.open(fn+".xyz", std::ofstream::out | std::ofstream::trunc);
        ofstr.close();

    }
    num_points = (int)(num_bp*density);
    point_dist = num_bp/num_points;

    min_peak_dist = PP_MIN_DETECT_DIST/chain->get_disc_len()/point_dist;
    loop_len      = PP_LOOP_LENGTH    /chain->get_disc_len()/point_dist;
    half_length   = loop_len/2;



//    sample_endpoints = arma::zeros<arma::ivec>(sample_points);

    std::cout << "Plasmid Endpoints: number of Points " << num_points << std::endl;
}
Dump_PlasmidEndpoints::~Dump_PlasmidEndpoints() {}


void Dump_PlasmidEndpoints::prod_dump() {
    call_count++;

    if (call_count%sample_step_dist==0 || sampling) {

        bool print_current_xyz=false;
        if (call_count%sample_step_dist==0) {
            sampling=true;
            sample_count=0;
            print_current_xyz=true;
        }
        sample_count++;

        // find reduced positions
        arma::mat pos(3,num_points);
        for (int i=0;i<num_points;i++) {
            pos.col(i) = chain->get_bp_pos()->col(i*point_dist);
        }

        arma::colvec writhe(num_points);
        arma::mat seg_pos(3,loop_len);
        for (int i=0;i<num_points;i++) {
            for (int j=0;j<loop_len;j++) {
                seg_pos.col(j) = pos.col(pmod(i-half_length+j,num_points));
            }
            writhe(i) = chain->cal_langowski_writhe_1a(seg_pos, false);
        }
        if (-arma::min(writhe)>arma::max(writhe)) {
            writhe = -writhe;
        }
        //writhe *= sgn(arma::mean(writhe));

        int smooth_box=PP_SMOOTHBOX/point_dist;
        arma::colvec smoothed = smooth(writhe,smooth_box);
        smoothed              = smooth(smoothed,smooth_box);

        /*
            Find first point outside of a peak
        */
//        sample_endpoints(sample_count-1) = find_num_endpoints(smoothed);;

//        vector<arma::colvec> endpoint_pos = find_endpoints_old(smoothed);
        std::vector<arma::colvec> endpoint_pos = find_endpoints(smoothed);

        if (sample_count==sample_points) {
            sampling=false;
//            ofstream file;
//            file.open(fn, ofstream::app);
//            file << find_endpoint_agreement(sample_endpoints) << "\n";
//            file.close();

            std::ofstream file;
            file.open(fn, std::ofstream::app);
            file << endpoint_pos.size();
            for (unsigned i=0;i<endpoint_pos.size();i++) {
                file << " " << endpoint_pos[i](0) << " " << endpoint_pos[i](1);
            }
            file << "\n";
            file.close();
        }

        if (call_count%sample_step_dist==0) {
            if (dump_writhe) {
                std::ofstream file;
                file.open(fn+".writhe", std::ofstream::app);
                for (int i=0;i<num_points;i++) {
                    file << smoothed(i) << ",";
                }
                file << smoothed(num_points-1) << "\n";
                file.close();
            }
            if (dump_xyz) {
                arma::colvec translate = arma::sum(*chain->get_bp_pos(),1)*1.0/num_bp;

                std::ofstream xyzfile;
                xyzfile.open(fn+".xyz", std::ofstream::app);
                xyzfile << num_points << " \n";
                xyzfile << "Atoms. Timestep: " << step << " \n";
                arma::vec spos;

//                bool endbead;
//                for (int i=0;i<num_points;i++) {
//                    spos = chain->get_bp_pos()->col(i*point_dist)-translate;
//                    endbead=false;
//                    for (int j=0;j<endpoint_pos.size();j++) {
//                        if (i==endpoint_pos[j](0) || pmod(i-1,num_points)==endpoint_pos[j](0) || pmod(i-2,num_points)==endpoint_pos[j](0) ||
//                                                     pmod(i+1,num_points)==endpoint_pos[j](0) || pmod(i+2,num_points)==endpoint_pos[j](0) ) {
//                            endbead=true;
//                            break;
//                        }
//                    }
//                    if (endbead) {
//                        xyzfile << "B " <<  spos(0) << " " << spos(1) << " " << spos(2) << " \n";
//                    }
//                    else {
//                        xyzfile << "A " <<  spos(0) << " " << spos(1) << " " << spos(2) << " \n";
//                    }
//                }

                int c = 0;
                for (int i=0;i<num_points;i++) {
                    spos = chain->get_bp_pos()->col(i*point_dist)-translate;
                    if (c%3==0) {
                        xyzfile << "A " <<  spos(0) << " " << spos(1) << " " << spos(2) << " \n";
                    }
                    if (c%3==1) {
                        xyzfile << "C " <<  spos(0) << " " << spos(1) << " " << spos(2) << " \n";
                    }
                    if (c%3==2) {
                        xyzfile << "N " <<  spos(0) << " " << spos(1) << " " << spos(2) << " \n";
                    }
                    c++;
                }
            }
        }
    }
}

//void Dump_PlasmidEndpoints::prod_dump() {
//
//    int half_length = loop_len/2;
//
//    // find reduced positions
//    arma::mat pos(3,num_points);
//    for (int i=0;i<num_points;i++) {
//        pos.col(i) = chain->get_bp_pos()->col(i*point_dist);
//    }
//
//    arma::colvec Writhe(num_points);
//    arma::mat seg_pos(3,loop_len);
//    for (int i=0;i<num_points;i++) {
//        for (int j=0;j<loop_len;j++) {
//            seg_pos.col(j) = pos.col(pmod(i-half_length+j,num_points));
//        }
//        Writhe(i) = chain->cal_langowski_writhe_1a(seg_pos, true);
//    }
//
//    ofstream file;
//    file.open(fn, ofstream::app);
//    for (int i=0;i<num_points;i++) {
//        file << Writhe(i) << ",";
//    }
//    file << Writhe(num_points-1) << "\n";
//    file.close();
//
//    ofstream xyzfile;
//    xyzfile.open(fn+".xyz", ofstream::app);
//    xyzfile << num_bp << " \n";
//    xyzfile << "Atoms. Timestep: " << step << " \n";
//    arma::vec spos;
//    for (int i=0;i<num_bp;i++) {
//        spos = chain->get_bp_pos()->col(i);
//        xyzfile << "A " <<  spos(0) << " " << spos(1) << " " << spos(2) << " \n";
//    }
//
//}




void Dump_PlasmidEndpoints::final_dump() {
//    prod_dump();
}

//////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////


arma::colvec Dump_PlasmidEndpoints::smooth(arma::colvec & writhe,int box_size) {
    arma::colvec smoothed(writhe.n_elem);

    double avg;
    for (unsigned i=0;i<writhe.n_elem;i++) {
        avg=0;
        for (int j=0;j<(2*box_size+1);j++) {
            avg+= writhe(pmod(i-box_size+j,writhe.n_elem));
        }
        smoothed(i) = avg/(2*box_size+1);
    }
    return smoothed;
}


std::vector<arma::colvec> Dump_PlasmidEndpoints::find_endpoints(arma::colvec & writhe) {

    std::vector<arma::colvec> end_points;
    std::vector<int> maxima;
    /*
        Find maximum as starting point
    */
    double max_val=0;
    int    max_id;
    for (int id=0;id<num_points;id++) {
        if (writhe(id) > max_val) {
            max_id  = id;
            max_val = writhe(id);
        }
    }
    if (max_val>PP_PEAK_THRESHOLD) {
        maxima.push_back(max_id);
    }
    else {
        return end_points;
    }
    /*
        Find all local maxima above the threshold
    */
    double local_min = max_val;
    int id;
    for (int i=1;i<num_points;i++) {
        id = pmod(max_id+i,num_points);
        if (writhe(id)>PP_PEAK_THRESHOLD) {
            /*
                Check for local maximum by comparing to the two neighboring values
            */
            if (writhe(pmod(id-1,num_points))<writhe(id) && writhe(pmod(id+1,num_points))<writhe(id)) {
                if (pmod(id-maxima.back(),num_points)<min_peak_dist) {
                    /*
                        If the two peaks are too close and the the new peak is larger the
                        last maximum is overwritten
                    */
                    if (writhe(id)>writhe(maxima.back())) {
                        maxima[maxima.size()-1]=id;
                        local_min = writhe(id);
                    }
                }
                else {
                    if (writhe(id)<=writhe(maxima.back())) {
                        if (writhe(id)-local_min > PP_LOCAL_MIN_DIFF) {
                            maxima.push_back(id);
                            local_min = writhe(id);
                        }
                    }
                    else {
                        if (writhe(maxima.back())-local_min > PP_LOCAL_MIN_DIFF) {
                            maxima.push_back(id);
                            local_min = writhe(id);
                        }
                        else {
                            maxima[maxima.size()-1]=id;
                            local_min = writhe(id);
                        }
                    }
                }
            }
        }
        if (writhe(id)<local_min) {
            local_min = writhe(id);
        }
    }
    /*
        Check if last maximum satisfies criteria
        In the case of a single maximum above the threshold this leads to the removal of that maximum.
        This is not a problem since this configuration should be discarded anyway.
    */
    if (pmod(max_id-maxima.back(),num_points)<min_peak_dist || writhe(maxima.back())-local_min < PP_LOCAL_MIN_DIFF) {
        maxima.erase(maxima.end() - 1);
    }

    for (unsigned i=0;i<maxima.size();i++) {
        end_points.push_back({double(maxima[i]),writhe(maxima[i])});
    }
    return end_points;
}


std::vector<arma::colvec> Dump_PlasmidEndpoints::find_endpoints_old(arma::colvec & smoothed) {
    int start = num_points + 1;
    int down_count=0;
    for (int i=0;i<num_points;i++) {
        if (smoothed(i) < PP_PEAK_THRESHOLD) {
            down_count++;
            if (down_count==PP_MIN_DETECT_DIST) {
                start = i;
                break;
            }
        }
        if (smoothed(i) > PP_PEAK_THRESHOLD) {
            down_count=0;
        }
    }

    // This should never occur!
    if (start > num_points) {
        std::cout << "start > num_points!" << std::endl;
        std::exit(0);
    }

    std::vector<arma::colvec> end_point_pos;

    int end_points=0;
    bool inpeak=false;
    bool possible_secondary_peak=false;
    double secondary_theshold=0;
    double last_peak=0;
    int    last_peak_pos=0;

    double val;

    for (int i=start;i<start+num_points;i++) {
        val = smoothed(pmod(i,num_points));
        if (smoothed(pmod(i,num_points))>PP_PEAK_THRESHOLD) {
            if (val>last_peak) {
                last_peak_pos = i;
                last_peak = val;
            }
//            if (!inpeak && down_count>=PP_MIN_DETECT_DIST) {
            if (!inpeak) {
                end_points++;
                inpeak=true;
//                down_count=0;
            }
        }
        if (val<PP_PEAK_THRESHOLD) {
            if (inpeak && !possible_secondary_peak) {
                end_point_pos.push_back({double(pmod(last_peak_pos,num_points)),last_peak});
            }
            inpeak=false;
            possible_secondary_peak=false;
            last_peak = 0;
//            down_count++;
        }

        if (inpeak && last_peak-val > PP_LOCAL_MIN_DIFF) {
            if (!possible_secondary_peak) {
                end_point_pos.push_back({double(pmod(last_peak_pos,num_points)),last_peak});
                secondary_theshold = val+PP_LOCAL_MIN_DIFF;
            }
            possible_secondary_peak=true;

            if (val+PP_LOCAL_MIN_DIFF < secondary_theshold) {
                secondary_theshold = val+PP_LOCAL_MIN_DIFF;
            }
        }
        if (possible_secondary_peak && val > secondary_theshold) {
            end_points++;
            last_peak_pos = i;
            last_peak = val;

            possible_secondary_peak=false;

            if (!inpeak) {
                std::cout << "Problem! " << std::endl;
                std::exit(0);
            }
        }

    }
    return end_point_pos;
}




int Dump_PlasmidEndpoints::find_endpoint_agreement(arma::ivec & num_endpoints) {
    int max_val = num_endpoints.max();
    arma::colvec counter(max_val+1);
    for (unsigned i=0;i<num_endpoints.n_elem;i++) {
        counter(num_endpoints(i))++;
    }
    counter=counter/num_endpoints.n_elem;
    for (int i=0;i<=max_val;i++) {
        if (counter(i) >= PP_MIN_AGREEMENT) {
            return i;
        }
    }
    return -1;
}









