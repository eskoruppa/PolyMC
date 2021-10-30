#include "Plasmid.h"


void plasmid_unwinding_sampling() {

//cout.precision(17);

    string interaction_file = "IDB/Plasmid1bp";
    string type             = "closed";

    long long int sweeps = 1e4;
    long long int steps  = 5e4;

    long long int equi_steps = 0;

    double T = 300;

    int num_bp    = 2600;
    string sequence = "a";

    double EV_rad     = 3.6;
    int num_EV_per_PB = 11;

    int dump_every = 2000;
//    string dump_dir = "dump/Branches/HT/NOBW_4kbp";
//    string dump_dir = "dump/Branches/HT/TEST";
    string dump_dir = "dump/Branches/test/2.6kbp_CS_#5";

    string dump_fn;

    double sigma = -0.05;

    double sigma_unwound    = -0.00;
    int steps_unwinding     = 5e3;
    int steps_winding       = 1e4;
    int step_equi           = 1e4;
    int n_split_winding     = 1;

    std::ostringstream strs;


    cout << "#############################################" << endl;
    cout << "#############################################" << endl;
    cout << "Initializing Simulation" << endl;
    PolyMC Run1;
    Run1.setup(type, interaction_file, num_bp, sequence, sigma, T);
    cout << "PolyMC - setup done" << endl;
    Run1.set_ExVol(EV_rad,num_EV_per_PB, true);
    cout << "PolyMC - ExVol initialized" << endl;
    strs.clear();
    strs.str("");
    strs << T;
    dump_fn = dump_dir + "_" + strs.str() + "K_";
    Run1.set_dumps(dump_fn,dump_every);
    cout << "PolyMC - Dumps initialized" << endl;

    ExVol* EV    = Run1.get_EV();
    Chain* chain = Run1.get_chain();

    Run1.unwinding_setup(EV, sigma_unwound, steps_unwinding, steps_winding, step_equi, n_split_winding);

    cout << "\n\n" << endl;
    cout << "#############################################" << endl;
    cout << "Equilibrating..." << endl;
    /*
        Equilibration
    */
    Run1.run(equi_steps,false, true);

    cout << "\n\n" << endl;
    cout << "#############################################" << endl;
    cout << "Commencing Production Run...\n" << endl;

    double loadtime;
    for (int sweep=0;sweep<sweeps;sweep++) {
        cout << "###########################" << endl;
        cout << "   sweep " << (sweep+1) << " \n\n";

//        cout << "Unwinding..." << endl;
//        Run1.unwinding_run();

//        cout << "Sampling..." << endl;
//        Run1.run(20,true, true);

//        cout << "###########################" << endl;
//        cout << "###########################" << endl;
//        cout << "###########################" << endl;
//        Run1.set_T(800);
//        cout << "T = 800..." << endl;
//        Run1.run(100000,true, true);
//
//        cout << "###########################" << endl;
//        cout << "###########################" << endl;
//        cout << "###########################" << endl;
//        Run1.set_T(700);
//        cout << "T = 300..." << endl;
//        Run1.run(50000,true, true);
//
//        cout << "###########################" << endl;
//        cout << "###########################" << endl;
//        cout << "###########################" << endl;
//        Run1.set_T(600);
//        cout << "T = 300..." << endl;
//        Run1.run(50000,true, true);
//
//        cout << "###########################" << endl;
//        cout << "###########################" << endl;
//        cout << "###########################" << endl;
//        Run1.set_T(500);
//        cout << "T = 300..." << endl;
//        Run1.run(50000,true, true);
//
//        cout << "###########################" << endl;
//        cout << "###########################" << endl;
//        cout << "###########################" << endl;
//        Run1.set_T(400);
//        cout << "T = 300..." << endl;
//        Run1.run(50000,true, true);


        cout << "###########################" << endl;
        cout << "###########################" << endl;
        cout << "###########################" << endl;
//        Run1.set_T(T);
//        cout << "T = " << T << "..." << endl;
        Run1.run(200000,true, true);

//        cout << "Sampling..." << endl;
//        Run1.run(steps,true, false);
    }
}



//void plasmid_unwinding_sampling() {
//
////cout.precision(17);
//
//    string interaction_file = "IDB/Plasmid6bp";
//    string type             = "closed";
//
//    double T = 300;
//
//    long long int sweeps = 1e5;
//    long long int steps  = 1e4;
//
//    long long int unwinding_steps     = 2e4;
//    long long int winding_steps       = 2e4;
//    long long int equilibration_steps = 1e5;
//    long long int sampling_steps      = 1e5;
//
//    double sigma                = -0.05;
//    double unwound_sigma        = -0.02;
//    double winding_increment    = -0.005;
//
//    vector<double> winding_sigmas;
//
//    double w_sigma=unwound_sigma+winding_increment;
//    while (std::abs(w_sigma)<std::abs(sigma)) {
//        winding_sigmas.push_back(w_sigma);
//        w_sigma+=winding_increment;
//    }
//
////    winding_sigmas.push_back(-0.005);
////    winding_sigmas.push_back(-0.010);
////    winding_sigmas.push_back(-0.015);
////    winding_sigmas.push_back(-0.020);
////    winding_sigmas.push_back(-0.025);
////    winding_sigmas.push_back(-0.030);
////    winding_sigmas.push_back(-0.035);
////    winding_sigmas.push_back(-0.040);
////    winding_sigmas.push_back(-0.045);
////    winding_sigmas.push_back(-0.050);
//
//
//    int num_bp    = 336;
//    string sequence = "a";
//
//    double EV_rad     = 3.51;
//    int num_EV_per_PB = 11;
//
//    string dump_dir = "dump/Branches/US/2.6k#3";
//    string dump_fn;
//
//    std::ostringstream strs;
//
//
//    cout << "#############################################" << endl;
//    cout << "#############################################" << endl;
//    cout << "Initializing Simulation" << endl;
//    PolyMC Run1;
//    Run1.setup(type, interaction_file, num_bp, sequence, sigma, T);
//    cout << "PolyMC - setup done" << endl;
//    Run1.set_ExVol(EV_rad,num_EV_per_PB, true);
//    cout << "PolyMC - ExVol initialized" << endl;
//    strs.clear();
//    strs.str("");
//    strs << T;
//    dump_fn = dump_dir + "_" + strs.str() + "K_";
//    Run1.set_dumps(dump_fn);
//    cout << "PolyMC - Dumps initialized" << endl;
//
//    ExVol* EV    = Run1.get_EV();
//    Chain* chain = Run1.get_chain();
//
//    double loadtime;
//    for (int sweep=0;sweep<sweeps;sweep++) {
//        cout << "###########################" << endl;
//        cout << "sweep " << (sweep+1) << endl;
//
//
//        /*
//            Winding
//        */
////        cout << "Winding..." << endl;
////        for (int w=0;w<winding_sigmas.size();w++) {
////            cout << "sigma set to " << winding_sigmas[w] << endl;
////            chain->set_sigma(winding_sigmas[w]);
////            EV->set_backup_conf(chain->get_bp_pos(),chain->get_triads());
////            loadtime = Run1.run(winding_steps,false,false);
////        }
////
////        cout << "Equilibrate..." << endl;
////        chain->set_sigma(sigma);
////        EV->set_backup_conf(chain->get_bp_pos(),chain->get_triads());
////        loadtime = Run1.run(equilibration_steps,true,false);
//
//        /*
//            Sampling
//        */
////        cout << "Sampling..." << endl;
//        loadtime = Run1.run(sampling_steps,true,false);
//
//        /*
//            Unwinding
//        */
////        cout << "Unwinding..." << endl;
////        cout << "sigma set to " << unwound_sigma << endl;
////        chain->set_sigma(unwound_sigma);
////        EV->set_backup_conf(chain->get_bp_pos(),chain->get_triads());
////        loadtime = Run1.run(unwinding_steps,false,false);
//
//    }
//}
















