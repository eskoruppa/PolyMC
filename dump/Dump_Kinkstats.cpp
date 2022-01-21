#include "Dump_Kinkstats.h"

Dump_Kinkstats::Dump_Kinkstats(Chain * ch, int N_dump, const std::string& filename, bool append)
: Dump(ch,N_dump,filename,append)
{
    if (!app) {
        std::ofstream ofstr;
        ofstr.open(fn, std::ofstream::out | std::ofstream::trunc);
        ofstr.close();
    }

    fn_kinkpos = fn+"_kinkpos";
    if (!app) {
        std::ofstream ofstr;
        ofstr.open(fn_kinkpos, std::ofstream::out | std::ofstream::trunc);
        ofstr.close();
    }

    num_omit_boundary = 50;
    first_considered  = num_omit_boundary;
    last_considered   = num_bp-num_omit_boundary-1;

    neglect_below_kink_dist = 50;


    tancor_maxdist = 200;
    tancor_sum   = arma::zeros(tancor_maxdist);
    elastic_sum  = arma::zeros(tancor_maxdist);
    tancor_count = arma::zeros(tancor_maxdist);

    fn_tancos = fn+"_tancos";
//    if (!app) {
//        std::ofstream ofstr;
//        ofstr.open(fn_tancos, std::ofstream::out | std::ofstream::trunc);
//        ofstr.close();
//    }

}

Dump_Kinkstats::~Dump_Kinkstats() {}

void Dump_Kinkstats::prod_dump() {

    int num_kink_left  = 0;
    int num_kink_right = 0;

    if (dump_fn_kinkpos) {
        kinkpos.clear();
    }


    std::vector<double> resp;
    for (int i=0;i<num_bps;i++) {
        resp = BPS[i]->get_status_diag();
        if (resp[0] == -1) {
            num_kink_left++;
        }
        if (resp[0] == +1) {
            num_kink_right++;
            if (dump_fn_kinkpos) {
                kinkpos.push_back(i);
            }
        }
    }

//    rkf_sum   += 1.*num_kink_right/num_bps;
//    rkf_count += 1;
//    std::cout << rkf_sum/rkf_count << std::endl;

    std::ofstream ofstr;
    ofstr.open(fn, std::ofstream::out | std::ofstream::app);
    ofstr  << 1.*num_kink_left/num_bps << " " << 1.*num_kink_right/num_bps << "\n";
    ofstr.close();

//    double ext            = arma::dot(chain->get_force_dir(),chain->get_bp_pos()->col(num_bp-1)-chain->get_bp_pos()->col(0));
//    double energy_elastic = chain->extract_true_energy();
//    double energy_work    = arma::dot(chain->get_beta_force_vec(),chain->get_bp_pos()->col(num_bp-1)-chain->get_bp_pos()->col(0))/chain->get_beta();
//
//    std::ofstream ofstr;
//    ofstr.open(fn, std::ofstream::out | std::ofstream::app);
//    ofstr  << num_bps << " " << num_kink_left << " " << num_kink_right << " " << ext << " " << energy_elastic << " " << energy_work << "\n";
//    ofstr.close();

    if (dump_fn_kinkpos) {

        ofstr.open(fn_kinkpos, std::ofstream::out | std::ofstream::app);
        for (unsigned i=0;i<kinkpos.size();i++) {
            if (i > 0) {
                ofstr << " ";
            }
            ofstr << kinkpos[i];
        }
        ofstr << "\n";
        ofstr.close();


        ///////////////////////////////////////////////////


        if (dump_fn_tancos) {

            arma::colvec tan,forcedir;
            forcedir = chain->get_force_dir();
            double tanfdot;
            double elastic_energy;
            for (unsigned i=first_considered;i<=last_considered;i++) {

                // find kink distance
                for (unsigned j=0;j<kinkpos.size();j++) {
                    unsigned dist,dist1,dist2;

                    dist = tancor_maxdist;
                    // break if there is a kink at pos i
                    if (j == 0) {
                        if (i < kinkpos[0]) {
                            dist = kinkpos[0] - i;
                        }
                    }
                    else if (j == kinkpos.size()-1) {
                        if (kinkpos[j] < i) {
                            dist = i - kinkpos[j];
                        }
                    }
                    else{
                        if (kinkpos[j] < i && i < kinkpos[j+1]) {
                            if (kinkpos[j+1] - kinkpos[j] < neglect_below_kink_dist) {
                                continue;
                            }

                            dist1 = i - kinkpos[j];
                            dist2 = kinkpos[j+1] - i;
                            if (dist1 < dist2) {
                                dist = dist1;
                            }
                            else {
                                dist = dist2;
                            }
                        }
                    }
                    if (dist < tancor_maxdist) {
                        dist--;
                        tan = chain->get_triads()->slice(i).col(2);
                        tanfdot = arma::dot(tan,forcedir);

                        (*chain->get_BPS())[i]->energy_extract_select();
                        elastic_energy = (*chain->get_BPS())[i]->energy_extract();

                        elastic_sum[dist]  += elastic_energy;
                        tancor_sum[dist]   += tanfdot;
                        tancor_count[dist] += 1;
                        break;
                    }
                }
            }
        }
    }
}

void Dump_Kinkstats::final_dump() {

	std::ofstream ofstr;
	if (app) ofstr.open(fn_tancos, std::ofstream::out | std::ofstream::app);
	else     ofstr.open(fn_tancos, std::ofstream::out | std::ofstream::trunc);

	double val;

    unsigned last=0;
    for (int n=0;n<tancor_sum.size();n++) {
        if (tancor_count[n]>0) {
            last++;
        }
        else {
            break;
        }
    }

    ofstr << last;
    for (int n=0;n<last;n++) {
        val = 1.*tancor_sum[n]/tancor_count[n];
        ofstr << " " << val;
    }
    ofstr << std::endl;

    ofstr << last;
    for (int n=0;n<last;n++) {
        val = 1.*elastic_sum[n]/tancor_count[n];
        ofstr << " " << val;
    }
    ofstr << std::endl;

    ofstr.close();
}




