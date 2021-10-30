#include "Plasmid.h"

void run_plasmid(int argc, const char **argv) {

    std::cout << std::endl;
    std::cout << "#############################################" << std::endl;
    std::cout << "######### Running Plasmid Protocol ##########" << std::endl;
    std::cout << "#############################################" << std::endl;
    std::cout << std::endl;
    std::string mode = "Plasmid";

    std::string IDB_fn = "IDB/Plasmid/TWLC10bp";

    long long int steps = 5e10;
    long long int equi  = 0;
    double T = 300;
    double sigma  = 0.00;
    int num_bp    = 200;
//    double f = 2.5;
//    arma::colvec fdir = {0,0,1};
    std::string dump_dir = "dump/Plasmid/";
    int  print  = 100000;
    bool Lk_check              = true;
    bool EV_active             = true;
    bool kill_if_link_voilated = false;
    bool check_energy_consistency = false;
    bool teardrop              = false;
    bool teardrop2d            = false;

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
    teardrop    = parse_flag("-teardrop",argc,argv);
    teardrop2d  = parse_flag("-teardrop2d",argc,argv);

    std::string seq="a";
    if (teardrop || teardrop2d) {
        seq = "b";
        for (int i=1;i<num_bp-1;i++) {
            seq += "a";
        }
        seq += "b";
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

    MCS_CSrot* rot      = new MCS_CSrot(&chain,2,larger(2,num_bp/2));
    MCS_Twist* twi      = new MCS_Twist(&chain);

    MCS_CSrot2d* rot2d    = new MCS_CSrot2d(&chain,2,larger(2,num_bp/4),{0,0,1});
    MCS_Slither2d* sli2d  = new MCS_Slither2d(&chain,2,larger(2,num_bp/8),{0,0,1});

    if (teardrop2d) {
//        MCSteps.push_back(rot2d);
        MCSteps.push_back(sli2d);
    }
    else {
        MCSteps.push_back(rot);
    }
    unsigned N_twist_moves = 1;
    for (unsigned i=0;i<N_twist_moves;i++) {
        MCSteps.push_back(twi);
    }



////// EXCLUDED VOLUMES /////
    double bp_per_EV  = 10;
    double EV_rad     = 0.34*(bp_per_EV+0.1);
    EV_rad    = parse_arg(EV_rad, "-EV"    ,argc,argv);
    if (EV_rad <= 0) {
        EV_active=false;
        Lk_check =false;
        EV_rad   = chain.get_disc_len();
    }
    ExVol EV(&chain,EV_rad);
    EV.set_repulsion_plane(false);

    if (EV_active) {
        for (unsigned i=0;i<MCSteps.size();i++) {
            MCSteps[i]->set_excluded_volume(&EV);
        }
    }
    for (int i=0;i<1e9;i++) {
        std::cout << "check repulsion plane option!" << std::endl;
    }


    std::vector<Dump*>   Dumps;

    /*
        Dump Writhe Map
    */
    int    WM_dump_every    = 0;
    double WM_seg_size       = 3.4;
    std::string WM_filename      = dump_dir+".wm";
    WM_dump_every  = parse_arg(WM_dump_every, "-WMn"  ,argc,argv);
    WM_seg_size    = parse_arg(WM_seg_size  , "-WMseg"  ,argc,argv);
    WM_filename    = parse_arg(WM_filename  , "-WMfn"  ,argc,argv);
    if (WM_dump_every>0) {
        Dump_WritheMap* DWM = new Dump_WritheMap(&chain, WM_dump_every, WM_filename, WM_seg_size, false);
        Dumps.push_back(DWM);
    }


    /*
        Dump Distance Map
    */

    int    DM_dump_every    = 0;
    double DM_density       = 0.5;
    std::string DM_filename      = dump_dir+".dm";
    DM_dump_every  = parse_arg(DM_dump_every, "-DMn"  ,argc,argv);
    DM_density     = parse_arg(DM_density   , "-DMseg"  ,argc,argv);
    DM_filename    = parse_arg(DM_filename  , "-DMfn"  ,argc,argv);
    if (DM_dump_every>0) {
        Dump_DistMap* DDM = new Dump_DistMap(&chain, DM_dump_every, DM_filename, DM_density, false);
        Dumps.push_back(DDM);
    }

    /*
        Dump State
    */
    int    State_dump_every = 0;
    std::string State_filename   = dump_dir+".state";

    State_dump_every        = parse_arg(State_dump_every, "-Stn"  ,argc,argv);
    State_filename          = parse_arg(State_filename  , "-Stfn" ,argc,argv);

    std::cout << State_filename << std::endl;
    if (State_dump_every>0) {
        Dump_State* DState = new Dump_State(&chain,mode, &EV, State_dump_every, State_filename ,false,false);
        Dumps.push_back(DState);
    }

    /*
        Dump Thetas
    */
    int    Thetas_dump_every = 0;
    std::string Thetas_filename   = dump_dir+".thetas";
    Thetas_dump_every        = parse_arg(Thetas_dump_every, "-Thn"  ,argc,argv);
    Thetas_filename          = parse_arg(Thetas_filename  , "-Thfn" ,argc,argv);

    if (Thetas_dump_every>0) {
        Dump_Thetas* DThetas = new Dump_Thetas(&chain, Thetas_dump_every, Thetas_filename ,false);
        Dumps.push_back(DThetas);
    }

    /*
        Dump Energy
    */
    int    Energy_dump_every = 0;
    std::string Energy_filename   = dump_dir+".en";
    Energy_dump_every        = parse_arg(Energy_dump_every, "-En"  ,argc,argv);
    Energy_filename          = parse_arg(Energy_filename  , "-Efn" ,argc,argv);

    if (Energy_dump_every>0) {
        Dump_Energy* DEnergy = new Dump_Energy(&chain, Energy_dump_every, Energy_filename ,false,false);
        Dumps.push_back(DEnergy);
    }

    /*
        Dump XYZ
    */
    int    XYZ_dump_every      = 0;
    std::string XYZ_dumpxyzfilename = dump_dir+".xyz";
    std::string XYZ_dumpxyzcenter   = "COM";
    std::string XYZ_rep             = "fl"; //"EV"
    XYZ_dump_every             = parse_arg(XYZ_dump_every,        "-XYZn"  ,argc,argv);
    XYZ_dumpxyzfilename        = parse_arg(XYZ_dumpxyzfilename  , "-XYZfn" ,argc,argv);
    if (XYZ_dump_every>0) {
        Dump_xyz* Dxyz = new Dump_xyz(&chain, XYZ_dump_every, XYZ_dumpxyzfilename ,false, XYZ_dumpxyzcenter,XYZ_rep);
        Dumps.push_back(Dxyz);
    }

    /*
        Dump format transformable to pdb
    */
    int    PDB_dump_every      = 0;
    std::string PDB_dumpxyzfilename = dump_dir+".triads";
    std::string PDB_dumpxyzcenter   = "COM";
    std::string PDB_rep             = "pdb"; //"EV"
    PDB_dump_every             = parse_arg(PDB_dump_every,        "-PDBn"  ,argc,argv);
    PDB_dumpxyzfilename        = parse_arg(PDB_dumpxyzfilename  , "-PDBfn" ,argc,argv);
    if (PDB_dump_every>0) {
        Dump_xyz* Dpdb = new Dump_xyz(&chain, PDB_dump_every, PDB_dumpxyzfilename ,false, PDB_dumpxyzcenter,PDB_rep);
        Dumps.push_back(Dpdb);
    }



    ///// DUMPS
//    Dump_Stiff              Dstiff(&chain, stiff_dump_every, stiff_filename,stiff_all);
//    Dump_EffTorsStiff       DETS(&chain, ETS_dump_every, ETS_filename, {0,0,1},ETS_writhe ,0.2*num_bp,0.8*num_bp,true);
//    Dump_ForceExtension     DFE(&chain, FE_dump_every, FE_filename,0.2*num_bp,0.8*num_bp,true);



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
            std::cout << "rot   acceptance rate = " << rot->acceptance_rate()   << std::endl;
            std::cout << "rot2d acceptance rate = " << rot2d->acceptance_rate()   << std::endl;
            std::cout << "sli2d acceptance rate = " << sli2d->acceptance_rate()   << std::endl;
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
