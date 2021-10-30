#include "Sim2d.h"

void run_teardrop2d(int argc, const char **argv) {

    std::cout << std::endl;
    std::cout << "#############################################" << std::endl;
    std::cout << "####### Running Teardrop2d Protocol #########" << std::endl;
    std::cout << "#############################################" << std::endl;
    std::cout << std::endl;
    std::string mode = "Teardrop2d";

    std::string IDB_fn = "IDB/Plasmid/TWLC10bp";

    long long int steps = 5e10;
    long long int equi  = 0;
    double T = 300;
    double sigma  = 0.00;
    int num_bp    = 200;
//    double f = 2.5;
//    arma::colvec fdir = {0,0,1};
    std::string dump_dir = "dump/Teardrop/";
    int  print  = 100000;
    bool Lk_check              = false;
    bool EV_active             = true;
    bool kill_if_link_voilated = false;
    bool check_energy_consistency = false;
    bool noflexhinge           = true;

    /*
        Parse Args
    */
    IDB_fn      = parse_arg(IDB_fn  , "-IDB"    ,argc,argv);
    steps       = parse_arg(steps   , "-steps"  ,argc,argv);
    equi        = parse_arg(equi    , "-equi"   ,argc,argv);
    T           = parse_arg(T       , "-T"      ,argc,argv);
    sigma       = parse_arg(sigma   , "-sigma"  ,argc,argv);
    num_bp      = parse_arg(num_bp  , "-nbp"    ,argc,argv);
    dump_dir    = parse_arg(dump_dir, "-dir"    ,argc,argv);
    noflexhinge = parse_flag("-noflexhinge",argc,argv);
//    teardrop2d  = parse_flag("-teardrop2d",argc,argv);

    int hingeelem = 12;

    int halfhing = 0.5*hingeelem;
    num_bp += halfhing*2;
    std::string hingeside = "";
    for (int i=0;i<halfhing;i++) {
        hingeside+="b";
    }

    std::string seq="a";
    if (!noflexhinge) {
        seq = hingeside;
        for (int i=halfhing;i<num_bp-halfhing;i++) {
            seq += "a";
        }
        seq += hingeside;
    }

    Chain chain(IDB_fn);
    chain.gen_circular(num_bp, sigma, seq);
    chain.fix_link();
    chain.set_T0_subtract(false);
    chain.set_T(T);

/*
    Monte Carlo Steps
*/

    std::vector<MCStep*> MCSteps;

    MCS_Twist* twi       = new MCS_Twist(&chain);
    MCS_CSrot2d* rot2d   = new MCS_CSrot2d(&chain,2,larger(2,num_bp/4),{0,0,1});
    MCS_CSrot_2dpot* rot2dpot = new MCS_CSrot_2dpot(&chain,2,larger(2,num_bp/4));
    MCS_Slither2d* sli2d = new MCS_Slither2d(&chain,2,larger(2,num_bp/4),{0,0,1});

//    MCSteps.push_back(rot2d);
    MCSteps.push_back(sli2d);
//    MCSteps.push_back(rot2dpot);

    int N_twist_moves = 0;
    for (int i=0;i<N_twist_moves;i++) {
        MCSteps.push_back(twi);
    }


////// EXCLUDED VOLUMES /////
    double bp_per_EV  = 10;
    double EV_rad     = 0.34*(bp_per_EV+0.1);
    EV_rad = 0;
    EV_rad    = parse_arg(EV_rad, "-EV"    ,argc,argv);
    if (EV_rad <= 0) {
        EV_active=false;
        Lk_check =false;
        EV_rad   = chain.get_disc_len();
    }
    ExVol EV(&chain,EV_rad);
    EV.set_repulsion_plane(true);


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
    std::cout << "Equilibration finished!" << std::endl;

    for (long long int step=0;step<=steps;step++) {

        if (step%print==0) {

            std::cout << "################" << std::endl;
            std::cout << "step " << step << std::endl;
            std::cout << "rot2d acceptance rate = " << rot2d->acceptance_rate()   << std::endl;
            std::cout << "sli2d acceptance rate = " << sli2d->acceptance_rate()   << std::endl;
            std::cout << "rotzp acceptance rate = " << rot2dpot->acceptance_rate()   << std::endl;
            std::cout << "twist acceptance rate = " << twi->acceptance_rate()   << std::endl;

            timer_finish = std::chrono::high_resolution_clock::now();
            timer_elapsed       = timer_finish - timer_start;
            std::cout << "elapsed time:   " << timer_elapsed.count() << " s\n";

            timer_start = std::chrono::high_resolution_clock::now();

            if (Lk_check) {
//                double Wr = chain.cal_langowski_writhe_1a();
                arma::mat writhemat;
                chain.langowski_writhe_elements(&writhemat);

                timer_finish = std::chrono::high_resolution_clock::now();
                timer_elapsed       = timer_finish - timer_start;
                std::cout << "WM timer:   " << timer_elapsed.count() << " s\n";
                timer_start = std::chrono::high_resolution_clock::now();

                double Wr = arma::accu(writhemat);
                double Tw = chain.cal_twist(0,num_bp);
                std::cout << "Wr = " << Wr << std::endl;
                std::cout << "Tw = " << Tw << std::endl;
                std::cout << "LK = " << Tw+Wr << " ( " << chain.get_dLK() << ")" << std::endl;

                if (std::abs(Tw+Wr-chain.get_dLK())>1) {
                    std::cout << "Linking Number not conserved" << std::endl;
                    if (kill_if_link_voilated) {
                        std::cout << "Terminating Simulation!" << std::endl;
                        std::exit(1);
                    }
                }

//                arma::mat writhemat;
//                chain.langowski_writhe_elements(&writhemat);
//                double Wr2 = arma::accu(writhemat);
//                std::cout << "LK2= " << Tw+Wr2 << " ( " << chain.get_dLK() << ")" << std::endl;
//                std::cout << "LK3= " << Tw+chain.gauss_writhe(1) << " ( " << chain.get_dLK() << ")" << std::endl;
            }

            if (check_energy_consistency) {
                if (!chain.check_energy_consistency()) {
                    std::cout << "Energy wrong!" << std::endl;
                }
                else {
                    std::cout << "Energy correct!" << std::endl;
                }
            }

        }
        if (step%100000==0) {
            /*
            Check numerical consistency of positions and triads.
            */
            if (!chain.config_consistent()) {
                std::cout << "chain positions inconsistent!" << std::endl;
                std::exit(0);
            }

            arma::mat  backup_pos    = *chain.get_bp_pos();
            arma::cube backup_triads = *chain.get_triads();

            chain.restore_consistency();
            if (chain.get_bp_pos()->has_nan()) {
                std::cout << "restore_consistency induces nan values in pos" << std::endl;
                chain.set_config(&backup_pos, &backup_triads, true);
//                std::exit(0);
            }

            if (EV_active) {
                if (EV.check_overlap()) {
                    std::cout << "Consistency check lead to overlap!" << std::endl;
                    std::exit(0);
                    EV.revert_to_backup();
                }
                else {
                    EV.set_current_as_backup();
                }
            }
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



void run_linear2d(int argc, const char **argv) {

    std::cout << std::endl;
    std::cout << "#############################################" << std::endl;
    std::cout << "######## Running Linear2d Protocol ##########" << std::endl;
    std::cout << "#############################################" << std::endl;
    std::cout << std::endl;
    std::string mode = "Linear2d";

    std::string IDB_fn = "IDB/Plasmid/TWLC10bp";

    long long int steps = 5e10;
    long long int equi  = 0;
    double T = 300;
    double sigma  = 0.00;
    int num_bp    = 200;
//    double f = 2.5;
//    arma::colvec fdir = {0,0,1};
    std::string dump_dir = "dump/Teardrop/";
    int  print  = 100000;
    bool Lk_check                   = false;
    bool EV_active                  = true;
    bool kill_if_link_voilated      = false;
    bool check_energy_consistency   = false;
    bool noflexhinge                = true;

    /*
        Parse Args
    */
    IDB_fn      = parse_arg(IDB_fn  , "-IDB"    ,argc,argv);
    steps       = parse_arg(steps   , "-steps"  ,argc,argv);
    equi        = parse_arg(equi    , "-equi"   ,argc,argv);
    T           = parse_arg(T       , "-T"      ,argc,argv);
    sigma       = parse_arg(sigma   , "-sigma"  ,argc,argv);
    num_bp      = parse_arg(num_bp  , "-nbp"    ,argc,argv);
    dump_dir    = parse_arg(dump_dir, "-dir"    ,argc,argv);
//    noflexhinge = parse_flag("-noflexhinge",argc,argv);

    std::string seq="a";

    Chain chain(IDB_fn);
    chain.gen_linear(num_bp,seq, sigma, {1,0,0} );
//    chain.fix_link();
    chain.set_T0_subtract(false);
    chain.set_T(T);

/*
    Monte Carlo Steps
*/

    std::vector<MCStep*> MCSteps;

    MCS_Twist* twi       = new MCS_Twist(&chain);
    MCS_CSrot2d* rot2d   = new MCS_CSrot2d(&chain,2,larger(2,num_bp/4),{0,0,1});
    MCS_Pivot2d* piv2d   = new MCS_Pivot2d(&chain,{0,0,1});

    MCS_Slither2d* sli2d = new MCS_Slither2d(&chain,2,larger(2,num_bp/4),{0,0,1},true);

//    MCSteps.push_back(rot2d);


    MCSteps.push_back(sli2d);
//    MCSteps.push_back(piv2d);

    int N_twist_moves = 0;
    for (int i=0;i<N_twist_moves;i++) {
        MCSteps.push_back(twi);
    }



////// EXCLUDED VOLUMES /////
    double bp_per_EV  = 10;
    double EV_rad     = 0.34*(bp_per_EV+0.1);
    EV_rad = 0;
    EV_rad    = parse_arg(EV_rad, "-EV"    ,argc,argv);
    if (EV_rad <= 0) {
        EV_active=false;
        Lk_check =false;
        EV_rad   = chain.get_disc_len();
    }
    ExVol EV(&chain,EV_rad);
    EV.set_repulsion_plane(true);


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
//            MCSteps[i]->MC();
            piv2d->MC();
        }
    }

    std::cout << "Equilibration finished!" << std::endl;

    for (long long int step=0;step<=steps;step++) {

        if (step%print==0) {

            std::cout << "################" << std::endl;
            std::cout << "step " << step << std::endl;
            std::cout << "rot2d acceptance rate = " << rot2d->acceptance_rate()   << std::endl;
            std::cout << "piv2d acceptance rate = " << piv2d->acceptance_rate()   << std::endl;
            std::cout << "sli2d acceptance rate = " << sli2d->acceptance_rate()   << std::endl;
            std::cout << "twist acceptance rate = " << twi->acceptance_rate()   << std::endl;

            timer_finish = std::chrono::high_resolution_clock::now();
            timer_elapsed       = timer_finish - timer_start;
            std::cout << "elapsed time:   " << timer_elapsed.count() << " s\n";

            timer_start = std::chrono::high_resolution_clock::now();

            if (check_energy_consistency) {
                if (!chain.check_energy_consistency()) {
                    std::cout << "Energy wrong!" << std::endl;
                }
                else {
                    std::cout << "Energy correct!" << std::endl;
                }
            }

        }

//        if (step%100000==0) {
//            piv2d->MC();
//        }
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



