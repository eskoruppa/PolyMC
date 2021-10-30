#include "open_elasticity.h"

void run_open(int argc, const char **argv) {

//std::cout.precision(17);

//    string IDB_fn = "IDB/Ceff_MS";
//    string IDB_fn = "IDB/Ceff_IB_Olson";
    std::string IDB_fn = "IDB/Ceff/Olson";
    std::string mode = "open";

    long long int steps = 1e7;
    long long int equi  = 0;

    double T = 300;

    double sigma  = 0.00;
    int num_bp    = 500;

    double f = 3.0;
    std::string dump_dir = "dump/open_elasticity/IB_Olson/IB";



    int print  = 100000;
    int cal_Lk = 1;
    bool Lk_first = false;

    bool EV_active = false;
    bool kill_if_link_voilated = false;


    /*
        Parse Args
    */

    IDB_fn      = parse_arg(IDB_fn  , "-IDB"    ,argc,argv);

    steps       = parse_arg(steps   , "-steps"  ,argc,argv);
    equi        = parse_arg(equi    , "-equi"   ,argc,argv);
    T           = parse_arg(T       , "-T"      ,argc,argv);
    sigma       = parse_arg(sigma   , "-sigma"  ,argc,argv);
    num_bp      = parse_arg(num_bp  , "-nbp"    ,argc,argv);
    f           = parse_arg(f       , "-f"      ,argc,argv);
    dump_dir    = parse_arg(dump_dir, "-dir"    ,argc,argv);




    Chain chain(IDB_fn);

    if (num_bp==500) {
//         chain.gen_linear(num_bp, "atttcccccagtgacttctggagacacgatgtacatcttctgatgcctagccttgatcactccgttatacctgtcccgtttcaatgcatcctaggcacgtaggaacatttttatacattaacctacagctaacaagaactctcgcatccacacgcgtttcttcccccggggaggcgaagcaagtttacgccggtcgcagagaactcacaagtagccgtgtccaacttcacctggcgccacgaccacagagtgggttgcagctagagtctacgctgtgaagaggataaggtccagagaaaactgacaaagatggatcacatatatttttattaacgacacatagaggttattatgtactgtagcgcaaacaagctaacccggcattcaacgctagtgcgacgaacgtggtatagttgtctcgcgcgatatgaaggccccgccggattggaggtgtccatccctgcgaactttccccatcgaacccgcacccgttactgccc", 0);
//         chain.gen_linear(num_bp, "cggtgacatatggccactgtcagtggtttaactccaactttttatattctatgaatctagacggcctatgtatcaggcgcaagtgaaggactgtggggattggcacgttaaccaggttttccttccttgaacactaatattttaggtctaactcaagcccactttctagacttacgcgacaaggagtgtgttctcgctgccataacggccaaacctgtatcgaaagatgcaattctctcctactccatacgtcttacgcaagaagctgtatgctgcggttagacgttattccattactataacttatacaatcatccatgaagtattcttcgtaaatcttaaccacgggtagcgtcactttgacttataacagcaattccgctgatttatcgaaaccgacgacaacgggctttctgtttatttgtaatacccctcaagacatgttatacgcaaacctaaccccgtagttataattgcttagagtcggactgtggggttag", 0);
         chain.gen_linear(num_bp, "acgagacggagccgtccggttgtgattttgccatgtacgtctcaccggcctatacgtattgccatactttcccagttgtcaagtacacggtcatcggatgagggaaaggtaactcaacttattttcgccgaggtcaaatatctccagggacggtttcgcagggctcgcatgagtttatcacagggtagcattatatattgtgcgaaggtagacctacccgtcgagagttaggtatctttcgggctgattggttttattcaactctagtgtctgcaaacacacgtcacgagttccttttgcgacgagacatggccgagaaaggtacttccattatgcagcgggtagacgcgctcactcttacccgagccagtactgccgtgctacggcatgcagactaaagaagctggttaagaggagttgttaagcgcaaatagcccgtactcaaaatataatgggcacaaagggcccaataacggaccgcactttagtgcggtataggc", 0);
    }
    else {
        chain.gen_linear(num_bp, "a", 0);
    }

    chain.set_force(f,{0,0,1});

    chain.set_T0_subtract(false);
    chain.set_T(T);


/*
    Monte Carlo Steps
*/
    std::vector<MCStep*> MCSteps;
    MCS_Pivot* piv = new MCS_Pivot(&chain);
    MCSteps.push_back(piv);

////// EXCLUDED VOLUMES /////
    double bp_per_EV  = 10;
    double EV_rad     = 0; // 0.34*(bp_per_EV+0.1);
    EV_rad    = parse_arg(EV_rad, "-EV"    ,argc,argv);
    if (EV_rad <= 0) {
        EV_active=false;
//        Lk_check =false;
        EV_rad   = chain.get_disc_len();
    }
    ExVol EV(&chain,EV_rad);
    EV.set_repulsion_plane(false);

    if (EV_active) {
        for (unsigned i=0;i<MCSteps.size();i++) {
            MCSteps[i]->set_excluded_volume(&EV);
        }
    }

    /*
        DUMPS
    */
    std::vector<Dump*>   Dumps = init_dump_cmdargs(argc,argv,&chain,mode,EV_rad);
    std::cout << "Dumps Initialized!" << std::endl;

    auto timer_start = std::chrono::high_resolution_clock::now();
    auto timer_finish = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> timer_elapsed;

    for (long long int step=0;step<equi;step++) {
        if (step%100000==0) {
            std::cout << "Equilibration-step " << step << std::endl;
        }
        for (unsigned i=0;i<MCSteps.size();i++) {
            MCSteps[i]->MC();
        }
    }
    std::cout << "Equilibrated" << std::endl;



    for (long long int step=0;step<=steps;step++) {

        if (step%print==0) {

            std::cout << "################" << std::endl;
            std::cout << "step " << step << std::endl;
            std::cout << "piv  acceptance rate = " << piv->acceptance_rate()   << std::endl;

            timer_finish = std::chrono::high_resolution_clock::now();
            timer_elapsed       = timer_finish - timer_start;
            std::cout << "elapsed time:   " << timer_elapsed.count() << " s\n";

            timer_start = std::chrono::high_resolution_clock::now();

            if (!chain.config_consistent()) {
                std::cout << "chain positions inconsistent!" << std::endl;
                std::exit(0);
            }

//            double Wr = chain.cal_langowski_writhe_1a(0.05);
//            double Tw = chain.cal_twist(0,num_bp);
//            std::cout << "Wr = " << Wr << std::endl;
//            std::cout << "Tw = " << Tw << std::endl;

        }

        for (unsigned i=0;i<MCSteps.size();i++) {
            MCSteps[i]->MC();
        }
        for (unsigned i=0;i<Dumps.size();i++) {
            Dumps[i]->dump();
        }
    }

    for (unsigned i=0;i<Dumps.size();i++) {
        Dumps[i]->final_dump();
    }
}
