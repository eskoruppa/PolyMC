#include "PolyMC.h"


PolyMC::PolyMC() :
step(0)
{

}
PolyMC::~PolyMC() {
    for (unsigned i=0;i<MCSteps.size();i++) {
        delete MCSteps[i];
    }
    for (unsigned i=0;i<Dumps.size();i++) {
        delete Dumps[i];
    }
}


void PolyMC::setup(string type, string interaction_file, int num_bp, string sequence, double sigma, double T) {
    this->type = type;

    Chain * new_chain = new Chain(interaction_file);
    if (type=="open") {
        new_chain->gen_linear(num_bp, sequence , sigma,{0,0,1});
    }
    if (type=="closed") {
        new_chain->gen_circular(num_bp, sigma, sequence);
        new_chain->fix_link();

//        MCS_CSrot * rot1 = new MCS_CSrot(new_chain,2,40);
//        MCS_CSrot * rot2 = new MCS_CSrot(new_chain,larger(2,num_bp/32),larger(2,num_bp/16));
//        MCS_CSrot * rot3 = new MCS_CSrot(new_chain,larger(2,num_bp/16),larger(2,num_bp/8));
//        MCS_CSrot * rot4 = new MCS_CSrot(new_chain,larger(2,num_bp/8),larger(2,num_bp/4));
//        MCS_CSrot * rot5 = new MCS_CSrot(new_chain,larger(2,num_bp/4),larger(2,num_bp/2));

        MCS_CSrot * rot1 = new MCS_CSrot(new_chain,2,larger(2,num_bp/2));
        MCS_CSrot * rot2 = new MCS_CSrot(new_chain,2,larger(2,num_bp/2));
        MCS_CSrot * rot3 = new MCS_CSrot(new_chain,2,larger(2,num_bp/2));
        MCS_CSrot * rot4 = new MCS_CSrot(new_chain,2,larger(2,num_bp/2));
        MCS_CSrot * rot5 = new MCS_CSrot(new_chain,2,larger(2,num_bp/2));

        MCSteps.push_back(rot1);
        MCSteps.push_back(rot2);
        MCSteps.push_back(rot3);
        MCSteps.push_back(rot4);
        MCSteps.push_back(rot5);

        /*
            Twist Moves
        */
        MCS_Twist * tw = new MCS_Twist(new_chain);

        int N_twist_moves = 40;
        for (int t=0;t<N_twist_moves;t++){
            MCSteps.push_back(tw);
        }


        /*
            BranchWinder Moves
        */
//        double seg_size         = 13.8;
//        double max_dist         = 12;
//        double min_contour_dist = 20;

//        MCS_BranchWinder * bw1 = new MCS_BranchWinder( new_chain, seg_size, max_dist, min_contour_dist );
//        MCSteps.push_back(bw1);

//        seg_size         = 27.4;
//        max_dist         = 14;
//        min_contour_dist = 40;
//
//        MCS_BranchWinder * bw2 = new MCS_BranchWinder( new_chain, seg_size, max_dist, min_contour_dist );
//        MCSteps.push_back(bw2);

//        seg_size         = 40;
//        max_dist         = 12;
//        min_contour_dist = 40;
//
//        MCS_BranchWinder * bw3 = new MCS_BranchWinder( new_chain, seg_size, max_dist, min_contour_dist );
//        MCSteps.push_back(bw3);

//        seg_size         = 7.4;
//        max_dist         = 12;
//        min_contour_dist = 20;
//
//        MCS_BranchWinder * bw4 = new MCS_BranchWinder( new_chain, seg_size, max_dist, min_contour_dist );
//        MCSteps.push_back(bw4);



//        double seg_size_min = 7;
//        double seg_size_max = 18;
//        MCS_Slither4 * sl4 = new MCS_Slither4(new_chain, seg_size_min, seg_size_max);
//        MCSteps.push_back(sl4);

//        int rep_size_min = 1;
//        int rep_size_max = 1;
//        int rep_dist_min = 10;
//        int rep_dist_max = 200;
//        bool translate_twist = true;
//        MCS_RepTherm * rep = new MCS_RepTherm(new_chain,rep_size_min,rep_size_max,rep_dist_min,rep_dist_max,translate_twist);
//        MCSteps.push_back(rep);

    }
    new_chain->set_T0_subtract(false);
    new_chain->set_T(T);

    chain = new_chain;
    setup_finalized = false;
}

void PolyMC::finalize_setup() {
    if (!setup_finalized) {
        step = 0;
        if (EV_active) {
            for (unsigned i=0;i<MCSteps.size();i++) {
                MCSteps[i]->set_excluded_volume(EV);
            }
        }
        setup_finalized = true;
    }
}


Chain* PolyMC::get_chain() {
    return chain;
}

ExVol* PolyMC::get_EV() {
    return EV;
}


void PolyMC::set_ExVol(double EV_rad, bool repulsion_plane) {
    EV_active = true;
    EV = new ExVol(chain,EV_rad);
    EV->set_repulsion_plane(repulsion_plane);
}


void PolyMC::set_dumps(string dir, int dump_every) {

    int    XYZ_dump_every      = dump_every/20;
    string XYZ_dumpxyzfilename = dir+"traj.xyz";
    string XYZ_dumpxyzcenter   = "COM";
    string XYZ_rep             = "dna"; //"EV"

    int    EN_dump_every  = 10;
    string EN_filename    = dir+"Energy";

    int    DPEP_dump_every         = dump_every;
//        string DPEP_filename           = "dump/Branches/plasmid_endpoints_4.6k_#14";
    string DPEP_filename           = dir+"plasmid_endpoints";
    double DPEP_density            = 0.05;
//    double DPEP_nm_len_loop        = 340; //1000
    int    DPEP_sample_step_dist   = 1;
    int    DPEP_sample_points      = 1;
    bool   DPEP_dump_writhe        = true;
    bool   DPEP_dump_xyz           = true;

    int    LINK_dump_every  = dump_every/100;
    string LINK_filename    = dir+"link";
    string LINK_options     = "quick";

    ///// DUMPS
    Dump_xyz*              Dxyz = new Dump_xyz(chain, XYZ_dump_every, XYZ_dumpxyzfilename ,false, XYZ_dumpxyzcenter,XYZ_rep);

    Dump_Energy*           DEN   = new Dump_Energy(chain, EN_dump_every, EN_filename,false);
    Dump_PlasmidEndpoints* DPEP  = new Dump_PlasmidEndpoints(chain, DPEP_dump_every, DPEP_filename, DPEP_density,DPEP_sample_step_dist,DPEP_sample_points,DPEP_dump_writhe,DPEP_dump_xyz);
    Dump_Linkingnumber*    DLINK = new Dump_Linkingnumber(chain, LINK_dump_every, LINK_filename,LINK_options,false);


    Dumps.push_back(DEN);

    Dumps.push_back(Dxyz);
    if (dump_every>0) {
        Dumps.push_back(DPEP);
//        Dumps.push_back(DLINK);
    }
}

/*
    MUTATORS
*/

void PolyMC::set_T(double T) {
    chain->set_T(T);
}







double PolyMC::run(long long int steps,bool dump, bool print) {

    finalize_setup();

    auto timer_start = std::chrono::high_resolution_clock::now();
    auto timer_finish = std::chrono::high_resolution_clock::now();
    chrono::duration<double> timer_elapsed;

    for (long long int current_step=0;current_step<=steps;current_step++) {
        step++;
        if (print) {
            if (step%1000==0) {
                cout << "################" << endl;
                cout << "step " << step << endl;
                timer_finish = std::chrono::high_resolution_clock::now();
                timer_elapsed       = timer_finish - timer_start;
                cout << "elapsed time:   " << timer_elapsed.count() << " s\n";
                timer_start = std::chrono::high_resolution_clock::now();

                double Wr = chain->cal_langowski_writhe_1a(0.5);
                double Tw = chain->cal_twist(0,chain->get_num_bps());
                double Lk = Wr+Tw;
                cout << "Tw = " << Tw << endl;
                cout << "Wr = " << Wr << endl;
                cout << "Lk = " << Lk << endl;

                int maxdecimals=int_log10(MCSteps.size());
                int decimals;
//                for (unsigned i=0;i<MCSteps.size();i++) {
//                    decimals = int_log10(i+1);
//                    cout << "#" << i+1 << " ";
//                    for (unsigned d=0;d<maxdecimals-decimals;d++) {
//                        cout << " ";
//                    }
//                    MCSteps[i]->print_acceptance_rate();
//                }

            }
        }
        if (step%cal_Lk_every==0) {
            if (!chain->check_link_conservation(0.5)) {
                if (kill_if_EV_voilated) {
                    cout << "Conservation of Linking Number violated!" << endl;
                    std::exit(0);
                }
            }
        }

        for (unsigned i=0;i<MCSteps.size();i++) {
            MCSteps[i]->MC();
        }

        if (dump) {
            for (unsigned i=0;i<Dumps.size();i++) {
                Dumps[i]->dump();
            }
        }

    }
    timer_finish = std::chrono::high_resolution_clock::now();
    timer_elapsed       = timer_finish - timer_start;
//    cout << "elapsed time:   " << timer_elapsed.count() << " s\n";
    return timer_elapsed.count();
}


//////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////// UNWINDING MOVE ////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////



void   PolyMC::unwinding_setup(ExVol * EV, double sigma_unwound, int steps_unwinding, int steps_winding, int step_equi, int n_split_winding) {
    this->sigma_unwound     = sigma_unwound;
    this->steps_unwinding   = steps_unwinding;
    this->steps_winding     = steps_winding;
    this->step_equi         = step_equi;
    this->n_split_winding   = n_split_winding;

    MCSU = new MCS_Unwinder(chain,sigma_unwound,steps_unwinding,steps_winding,step_equi,n_split_winding );

    MCSU->set_excluded_volume(EV);
}


double PolyMC::unwinding_run() {

    auto timer_start = std::chrono::high_resolution_clock::now();
    auto timer_finish = std::chrono::high_resolution_clock::now();
    chrono::duration<double> timer_elapsed;

    if (MCSU->MC()) {
        cout << "Unwinding Move Accepted!" << endl;
    }
    else {
        cout << "Unwinding Move Rejected!" << endl;
    }

    timer_finish = std::chrono::high_resolution_clock::now();
    timer_elapsed       = timer_finish - timer_start;
//    cout << "elapsed time:   " << timer_elapsed.count() << " s\n";
    return timer_elapsed.count();
}








