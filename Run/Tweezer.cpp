#include "Tweezer.h"

/*

*/

void run_tweezer(int argc, const char **argv) {

    std::cout << std::endl;
    std::cout << "######################################" << std::endl;
    std::cout << "###### Running Tweezer Protocol ######" << std::endl;
    std::cout << "######################################" << std::endl;
    std::cout << std::endl;
    std::string mode = "Tweezer";


    std::string IDB_fn = "IDB/Tweezer/TWLC10bp";
//    string IDB_fn = "";

    long long int steps = 0;
    long long int equi  = 0;
    int num_bp          = 0;

    double T      = 300;
    double sigma  = 0.00;
    double f      = 0;
    double EV_rad = 0;
    arma::colvec fdir = {0,0,1};
    std::string dump_dir = "Backdoorsluts9";

    int  print_every  = 100000;
    bool EV_active = true;
    bool kill_if_link_voilated   = false;
    bool check_link_conservation = true;
    int  consistency_check_every = 100000;
    bool copy_input = false;


    std::string inputfn = "";
    inputfn  = parse_arg(inputfn, "-in",argc,argv);

    if (inputfn!="") {
        InputRead input(inputfn);
        input.testfilefound();

        IDB_fn      = input.get_single_val<std::string>   ("IDB");
        steps       = input.get_single_val<long long int> ("steps");
        equi        = input.get_single_val<long long int> ("equi",0);
        T           = input.get_single_val<double>        ("T",T);
        num_bp      = input.get_single_val<int>           ("num_bp");
        sigma       = input.get_single_val<double>        ("sigma",0);
        f           = input.get_single_val<double>        ("f", 0);
        EV_rad      = input.get_single_val<double>        ("EV",0);
        dump_dir    = input.get_single_val<std::string>   ("dump_dir");
        print_every = input.get_single_val<int>           ("print_every");
        copy_input  = input.get_single_val<bool>          ("copy_input",false);
    }

    /*
        Parse Args
    */
    IDB_fn      = parse_arg(IDB_fn  , "-IDB"    ,argc,argv);
    steps       = parse_arg(steps   , "-steps"  ,argc,argv);
    equi        = parse_arg(equi    , "-equi"   ,argc,argv);
    num_bp      = parse_arg(num_bp  , "-nbp"    ,argc,argv);
    T           = parse_arg(T       , "-T"      ,argc,argv);
    sigma       = parse_arg(sigma   , "-sigma"  ,argc,argv);
    f           = parse_arg(f       , "-f"      ,argc,argv);
    dump_dir    = parse_arg(dump_dir, "-dir"    ,argc,argv);



    /*
        Copy IDB and Input file to dump dir
    */
    if (copy_input) {
        std::string copy_inputfn = dump_dir+".in";
        std::string copy_IDBfn   = dump_dir+".idb";

        std::ifstream  src_in(inputfn,      std::ios::binary);
        std::ofstream  dst_in(copy_inputfn, std::ios::binary);
        dst_in << src_in.rdbuf();

        std::ifstream  src_idb(IDB_fn,     std::ios::binary);
        std::ofstream  dst_idb(copy_IDBfn, std::ios::binary);
        dst_idb << src_idb.rdbuf();
    }



    Chain chain(IDB_fn);
    chain.gen_linear(num_bp, "a", sigma,fdir);
//    chain.gen_linear(num_bp, "a", 0.005,fdir);

    chain.fix_termini_orientation();
    chain.fix_termini_position();
    chain.set_force(f,fdir);

    chain.set_T0_subtract(false);
    chain.set_T(T);

    arma::cube backup_triads = *chain.get_triads();
    arma::mat  backup_bp_pos = *chain.get_bp_pos();


    chain.impose_torsional_trap(chain.get_dLK());


//    chain.set_dLK_fix(chain.get_dLK());
//    chain.set_measure_torque(100);
//    chain.set_torsional_trap(chain.get_dLK());


    /*
        Monte Carlo Steps
    */

    std::vector<MCStep*> MCSteps;

//    MCS_Pivot* piv        = new MCS_Pivot(&chain);
//    MCS_CSrot* rot1       = new MCS_CSrot(&chain,larger(2,num_bp/4),larger(2,num_bp/2));
    MCS_CSrot* rot2       = new MCS_CSrot(&chain,2,larger(2,num_bp/4));
//    MCS_CStrans* trans    = new MCS_CStrans(&chain,2,num_bp/2,{0,0,1});
//    MCS_CStrans* trans    = new MCS_CStrans(&chain,2,num_bp/2,{0,0,0});
//    MCS_PivCon* pivcon    = new MCS_PivCon(&chain,2,num_bp/4);
    MCS_PivCon* pivcon    = new MCS_PivCon(&chain,2,num_bp/10);
//    MCS_Twist* twi        = new MCS_Twist(&chain);
    MCS_TailTwist* tailtw = new MCS_TailTwist(&chain);

//    MCS_torquePiv* torpiv = new MCS_torquePiv(&chain);


//    MCSteps.push_back(piv);

    MCSteps.push_back(rot2);
    MCSteps.push_back(pivcon);
    MCSteps.push_back(tailtw);

//    MCSteps.push_back(torpiv);

/*
    MCSteps.push_back(rot1);
    MCSteps.push_back(trans);
*/

//    int N_twist_moves = 0;
//    for (int i=0;i<N_twist_moves;i++) {
//        MCSteps.push_back(twi);
//    }


////// EXCLUDED VOLUMES /////
    EV_rad    = parse_arg(EV_rad, "-EV"    ,argc,argv);
    if (EV_rad <= 0) {
        EV_active=false;
        EV_rad   = chain.get_disc_len();
    }
    ExVol EV(&chain,EV_rad);
    EV.set_repulsion_plane(EV_active);

    if (EV_active) {
        for (unsigned i=0;i<MCSteps.size();i++) {
            MCSteps[i]->set_excluded_volume(&EV);
        }
    }

    /*
        DUMPS
    */
    std::vector<Dump*>   Dumps = init_dump_cmdargs(argc,argv,&chain,mode,EV_rad,inputfn);
    std::cout << "Dumps Initialized!" << std::endl;

    /*
        Equilibration
    */
    std::cout << "Starting Equilibration" << std::endl;
    for (long long int step=0;step<equi;step++) {
        if (step%100000==0) {
            std::cout << "Equilibration-step " << step << std::endl;
        }



        for (unsigned i=0;i<MCSteps.size();i++) {
            MCSteps[i]->MC();
        }
    }
    std::cout << "Equilibrated" << std::endl;
    std::cout << "starting simulations" << std::endl;

//    double tw_sum          = 0;
//    double wr_sum          = 0;
//    long long int tw_count = 0;


/*
    Slowly introduce sigma
*/

//    double dLK_target = chain.sigma2dLk(sigma);
//    chain.set_torstrap_dLK_fix(1);
//
//    long long int add_dlk_steps=0;
//    double dLK_fix = chain.get_torstrap_dLK_fix();
//    while (dLK_fix<dLK_target) {
//        add_dlk_steps++;
//        for (unsigned i=0;i<MCSteps.size();i++) {
//            MCSteps[i]->MC();
//        }
//
//        for (unsigned i=0;i<Dumps.size();i++) {
//            Dumps[i]->dump();
//        }
//        if (add_dlk_steps%50000==0) {
//
//            std::cout << "################" << std::endl;
//            std::cout << "step " << add_dlk_steps << std::endl;
//            std::cout << "fix  = " << chain.get_torstrap_dLK_fix() << std::endl;
//            std::cout << "LK   = " << chain.get_dLK() << " ( " <<  dLK_target <<  " )" << std::endl;
//
////            double Wr = chain.cal_langowski_writhe_1a(1);
////            double Tw = chain.cal_twist(0,num_bp);
////            std::cout << "Wr = " << Wr << std::endl;
//////                std::cout << "Wr2= " << chain.cal_langowski_writhe_1a(*chain.get_bp_pos(),false) << std::endl;
////            std::cout << "Tw = " << Tw << std::endl;
////            std::cout << "LK = " << Tw+Wr << " ( " << chain.get_dLK() << " / " <<  dLK_target <<  " )" << std::endl;
//        }
//        if (add_dlk_steps%50000==0) {
//            dLK_fix++;
//            chain.set_torstrap_dLK_fix(dLK_fix);
//        }
//
//    }


    /*
        Main Loop
    */

    auto timer_start = std::chrono::high_resolution_clock::now();
    auto timer_finish = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> timer_elapsed;

    double lk  = 0;
    double lk1 = 0;
    double lk2 = 0;
    long long lk_count = 0;

    double maxdiff     = 0;
    double maxdiff_add = 0;

    for (long long int step=1;step<=steps;step++) {

        if (step%print_every==0) {

            std::cout << "################" << std::endl;
            std::cout << "step " << step << std::endl;
//            for (int i=0;i<MCSteps.size();i++) {
//                MCSteps[i]->print_acceptance_rate();
//            }

//            std::cout << "piv    acceptance rate = " << piv->acceptance_rate()   << std::endl;
            std::cout << "pivco  acceptance rate = " << pivcon->acceptance_rate()   << std::endl;
//            std::cout << "tqpiv acceptance rate  = " << torpiv->acceptance_rate()   << std::endl;
//            std::cout << "rot1  acceptance rate = " << rot1->acceptance_rate()   << std::endl;
            std::cout << "rot2   acceptance rate = " << rot2->acceptance_rate()   << std::endl;
//            std::cout << "trans acceptance rate = " << trans->acceptance_rate()   << std::endl;
//            std::cout << "twist acceptance rate = " << twi->acceptance_rate()   << std::endl;
            std::cout << "tailtw acceptance rate = " << tailtw->acceptance_rate()   << std::endl;


//            std::cout << "Mean Torque = " << chain.get_torque_measure_accu(false)/(2*M_PI) << std::endl;
//            std::cout << "Current dLK = " << chain.get_dLK() << std::endl;
//            std::cout << "dLK fix     = " << chain.get_torstrap_dLK_fix() << std::endl;
//            std::cout << "dLK aim     = " << chain.get_torstrap_dLK_aim() << std::endl;
//            std::cout << "mean LK     = " << lk1/lk_count << std::endl;
//            std::cout << "len         = " << chain.get_contour_len() << std::endl;
//            std::cout << "Ceff        = " << chain.get_contour_len()/(4*M_PI*M_PI*(lk2/lk_count-lk1/lk_count*lk1/lk_count)) << std::endl;
//
//
            timer_finish = std::chrono::high_resolution_clock::now();
            timer_elapsed       = timer_finish - timer_start;
            std::cout << "elapsed time:   " << timer_elapsed.count() << " s\n";

            timer_start = std::chrono::high_resolution_clock::now();


            double Wr = chain.cal_langowski_writhe_1a(1);
            double Tw = chain.cal_twist(0,num_bp);
            std::cout << "Wr = " << Wr << std::endl;
//                std::cout << "Wr2= " << chain.cal_langowski_writhe_1a(*chain.get_bp_pos(),false) << std::endl;
            std::cout << "Tw = " << Tw << std::endl;
            std::cout << "LK = " << Tw+Wr << " ( " << chain.get_dLK() << ")" << std::endl;

//            double diff = chain.get_dLK() - Tw-Wr;
//            std::cout << "diff    = " << diff  << std::endl;
//
//            int add = 2;
//            arma::colvec zdir = {0,0,1};
//            double closfac    = 10;
//
//            arma::mat addpos = arma::zeros(3,num_bp+2*add);
////            addpos.submat( add, add, num_bp-1+add, num_bp-1+add) = *chain.get_bp_pos();
//
//            for (int i=0;i<num_bp;i++) {
//                addpos.col(i+add) = (*chain.get_bp_pos()).col(i);
//            }
//
//            for (int i=add-1;i>=0;i--) {
//                addpos.col(i) = addpos.col(i+1)-zdir*closfac;
//            }
//            for (int i=num_bp+add;i<num_bp+2*add;i++) {
//                addpos.col(i) = addpos.col(i-1)+zdir*closfac;
//            }
//
//            double Wr_add = chain.cal_langowski_writhe_1a(addpos,false);
//            std::cout << "LK_add   = " << Tw+Wr_add << " ( " << chain.get_dLK() << ")" << std::endl;
//            double diff_add = chain.get_dLK() - Tw-Wr_add;
//            std::cout << "diff_add = " << diff_add  << std::endl;
//
//            if (std::abs(diff) > maxdiff) {
//                maxdiff = std::abs(diff);
//            }
//            std::cout << "maxdiff     = " << maxdiff  << std::endl;
//            if (std::abs(diff_add) > maxdiff_add) {
//                maxdiff_add = std::abs(diff_add);
//            }
//            std::cout << "maxdiff_add = " << maxdiff_add  << std::endl;


//            for (int i=0;i<num_bp+2*add;i++) {
//                if (arma::norm(addpos.col(i))==0) {
//                    std::cout << i << std::endl;
//                    std::cout << addpos.col(i).t();
//                    if (i!=add ) {
//                        std::std::cout << "null col found" << std::endl;
//                        std::exit(0);
//                    }
//
//                }
//            }



        }

        /*
        if (step%consistency_check_every==0) {

            if (check_link_conservation) {
                double Wr = chain.cal_langowski_writhe_1a(1);
                double Tw = chain.cal_twist(0,num_bp);
                std::cout << "Wr = " << Wr << std::endl;
//                std::cout << "Wr2= " << chain.cal_langowski_writhe_1a(*chain.get_bp_pos(),false) << std::endl;
                std::cout << "Tw = " << Tw << std::endl;
                std::cout << "LK = " << Tw+Wr << " ( " << chain.get_link() << ")" << std::endl;

                int add = 2;
                arma::colvec zdir = {0,0,1};
                double closfac    = 10;

                arma::mat addpos = arma::zeros(3,num_bp+2*add);
                addpos.submat( add, add, num_bp-1+add, num_bp-1+add) = *chain.get_bp_pos();

                arma::colvec zdir = {0,0,1};
                double closfac    = 10;

                arma::mat cpos = *chain.get_bp_pos();
                for (int i=add-1;i>=0;i--) {
                    cpos.col(i) = cpos.col(i+1)-zdir*closfac;
                }
                for (int i=num_bp+add;i<num_bp+2*add;i++) {
                    cpos.col(i) = cpos.col(i-1)+zdir*closfac;
                }


                if (std::abs(chain.get_link()-Tw-Wr)> 1) {
                    if (kill_if_link_voilated) {
                        std::cout << "linking number violated!" << std::endl;
                        std::exit(0);
                    }
                    else {
                        std::cout << "Linking Number violated... reverting to backup configuration." << std::endl;
                        chain.set_config(&backup_bp_pos, &backup_triads, false, false);
                        EV.set_current_as_backup();

                        string lkviolatedfn = "dump/test/test.linkviolated";
                        lkviolatedfn = parse_arg(dump_dir, "-dir"    ,argc,argv)+".linkviolated";
                        ofstream ofstr;
                        ofstr.open(lkviolatedfn, ofstream::out | ofstream::app);
                        ofstr << step << " " << chain.get_link() << " " << Wr << " " << Tw << std::endl;
                        ofstr.close();
                    }
                }
                else {
                    backup_triads = *chain.get_triads();
                    backup_bp_pos = *chain.get_bp_pos();
                }
            }

            bool consistent = chain.check_energy_consistency();
            if (!consistent) {
                std::cout << "Energy inconsistency!" << std::endl;
//                std::exit(0);
            }

            if (!chain.config_consistent()) {
                std::cout << "chain positions inconsistent!" << std::endl;
                chain.restore_consistency();
                if (EV_active) {
                    if (EV.check_overlap()) {
                        std::cout << "Consistency check lead to overlap!" << std::endl;
                        EV.revert_to_backup();
                    }
                    else {
                        EV.set_current_as_backup();
                    }
                }
            }


            double zlast = chain.get_bp_pos()->col(num_bp-1)(2);
            for (int i=0;i<num_bp-1;i++) {
                if (chain.get_bp_pos()->col(i)(2) >= zlast) {
                    std::cout << "Repulsion Plane violated!" << std::endl;
                    string lkviolatedfn = "dump/test/test.linkviolated";
                    lkviolatedfn = parse_arg(dump_dir, "-dir"    ,argc,argv)+".linkviolated";
                    ofstream ofstr;
                    ofstr.open(lkviolatedfn, ofstream::out | ofstream::app);
                    ofstr << step << " Repulsion Plane Violated" << std::endl;
                    ofstr.close();
                }
            }

        }
        */




        /*
        if (step%1000==0) {
            int from = 0.25*num_bp;
            int to   = 0.75*num_bp;

            arma::colvec zdir = {0,0,1};
            double closfac    = 10;

//            arma::mat wrelem;
//            chain.langowski_writhe_elements(&wrelem);
//            double lwr = arma::accu(wrelem.submat( from, from, to, to ));

            arma::mat cpos = *chain.get_bp_pos();
            for (int i=from-1;i>=0;i--) {
                cpos.col(i) = cpos.col(i+1)-zdir*closfac;
            }
            for (int i=to+1;i<num_bp;i++) {
                cpos.col(i) = cpos.col(i-1)+zdir*closfac;
            }


            double lwr = chain.cal_langowski_writhe_1a(cpos,false);
            double fwr = chain.cal_fuller_writhe(from, to, zdir);

            string wrcomp = "dump/test/writhecomparision";

            ofstream ofstr;
            ofstr.open(wrcomp, ofstream::out | ofstream::app);
            ofstr << lwr << " " << fwr << std::endl;
            ofstr.close();
        }
        */







        for (unsigned i=0;i<MCSteps.size();i++) {
            MCSteps[i]->MC();
        }
        for (unsigned i=0;i<Dumps.size();i++) {
            Dumps[i]->dump();
        }


        lk  = chain.get_dLK();
        lk1 += lk;
        lk2 += lk*lk;
        lk_count++;
    }

    for (unsigned i=0;i<Dumps.size();i++) {
        Dumps[i]->final_dump();
    }

}
