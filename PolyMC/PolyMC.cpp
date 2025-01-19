#include "../PolyMC.h"


PolyMC::PolyMC(const std::vector<std::string> & argv, bool produce_output) :
input_given(false),nooutput(!produce_output),geninfile()
{
    std::cout << std::endl;
    std::cout << "######################################" << std::endl;
    std::cout << "######## Initializing PolyMC #########" << std::endl;
    std::cout << "######################################" << std::endl;
    std::cout << std::endl;
    init_geninfile();

//    this->argc = argc;
    this->argv = argv;

    mode = parse_str_atpos(1,argv);

    inputfn  = "";
    inputfn  = parse_arg(inputfn, "-in",argv);
    InputRead * inputread = new InputRead(inputfn);
    input = inputread;

    if (inputfn!= "") {
        if (!input->filefound()) {
            std::cout << "Error: Input file '" << inputfn << "' not found!" << std::endl;
            std::exit(0);
        }

        if (input->contains_single("mode")) {
            mode = input->get_single_val("mode",mode);
        }
        input_given = true;
    }

    std::cout << " mode             = " << mode << std::endl;
    geninfile.add_entry(GENINFILE_SIMSETUP,"mode",mode);

    init_general();
    init_seed();
    seq = read_seq(num_bp);

    #ifdef INCLUDE_OPEN_MODE
    if (mode == POLYMC_MODE_OPEN) {
        init_open();
    }
    #endif
    #ifdef INCLUDE_TWEEZER_MODE
    if (mode == POLYMC_MODE_TWEEZER) {
        init_tweezer();
    }
    #endif
    #ifdef INCLUDE_PLECTONEME_MODE
    if (mode == POLYMC_MODE_PLEC) {
        init_plec();
    }
    #endif
    #ifdef INCLUDE_UMBRELLAPLASMID_MODE
    if (mode == POLYMC_MODE_UMBRELLAPLASMID) {
        init_umbrellaplasmid();
    }
    #endif
    #ifdef INCLUDE_PLASMID_MODE
    if (mode == POLYMC_MODE_PLASMID) {
        init_plasmid();
    }
    #endif
    #ifdef INCLUDE_2D_MODES
    if (mode == POLYMC_MODE_LIN2D) {
        init_lin2d();
    }
    if (mode == POLYMC_MODE_CLOSED2D) {
        init_closed2d();
    }
    #endif


/*----------------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------------*/


//    Pair * p      = new Pair();
//    Pair_LJ * plj = new Pair_LJ();
//
//    p->test();
//    plj->test();

//    std::exit(0);

//    strvec2doublevec

//    Unbound * unb = new Unbound(chain);
//    unb->read_inputfile("UnboundInteractions");
//    unb->init_monomer_type_assignments_from_file("SetAtoms");
//
//    std::exit(0);

/*----------------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------------*/


    if (!chain_initialized) {
        if (mode == "") {
            std::cout << "No simulation mode specified!" << std::endl;
        }
        else {
            std::cout << "Error: Unknown simulation mode '" << mode << "'!" << std::endl;
        }
        std::cout << "Valid modes are:" << std::endl;
        std::cout << " - " << POLYMC_MODE_TWEEZER           << std::endl;
        std::cout << " - " << POLYMC_MODE_PLEC              << std::endl;
        std::cout << " - " << POLYMC_MODE_PLASMID           << std::endl;
        std::cout << " - " << POLYMC_MODE_UMBRELLAPLASMID   << std::endl;
        std::cout << " - " << POLYMC_MODE_OPEN              << std::endl;
        std::cout << " - " << POLYMC_MODE_LIN2D             << std::endl;
        std::cout << " - " << POLYMC_MODE_CLOSED2D          << std::endl;
        std::exit(0);
    }

    /*
        Initialialze additional interactions and constraints
    */

    init_ExVol();
    init_ElStat(ElStat_fn);
    init_constraints();
    init_constraints(constraint_fn, true);
    init_GenForce(GenForce_fn);
    init_pair_interactions();
    init_stateswitch_moves();

    if (!nooutput) {
        init_dump();
        /*
            Copy IDB and Input file to dump dir
        */
        if (copy_input) {
            copy_input_files();
        }
    }

    if (slow_wind) {
        if (mode == POLYMC_MODE_TWEEZER) {
            tweezer_slow_wind(0, chain->get_dLK(), slow_wind_dLK_step, slow_wind_steps_per_dLK_step);
        }
        else if (mode == POLYMC_MODE_PLEC) {
            plec_slow_wind   (0, chain->get_dLK(), slow_wind_dLK_step, slow_wind_steps_per_dLK_step);
        }
        else {
            generic_slow_wind(0, chain->get_dLK(), slow_wind_dLK_step, slow_wind_steps_per_dLK_step);
        }
    }

//    (*chain->get_BPS())[100]->deactivate_twist_energy();

    if (mode == POLYMC_MODE_LIN2D && use_slither2d) {
        lin2d_slither_equi_with_pivot();
    }
    
}


PolyMC::~PolyMC() {
    for (unsigned i=MCSteps.size()-1;i>=1;i--) {
        for (unsigned j=0;j<i;j++) {
            if (MCSteps[i]==MCSteps[j]) {
                MCSteps.erase(MCSteps.begin() + j);
                break;
            }
        }
    }
    for (unsigned i=0;i<MCSteps.size();i++) {
        delete MCSteps[i];
    }
    MCSteps.clear();
    for (unsigned i=0;i<Dumps.size();i++) {
        delete Dumps[i];
    }
    Dumps.clear();
    if (chain_initialized) {
        delete chain;
    }
    delete input;
    if (EV_active) {
        delete EV;
    }
    if (ES_active) {
        delete ES;
        delete ES_Potential;
    }
}


bool PolyMC::simple_execute() {
    if (equi>0) {
        run(equi,MCSteps,"Equilibration",false);
    }
    ///////////////////////////////////////////////
    // Main loop
    run(steps,MCSteps,"Production Run");
    ///////////////////////////////////////////////
    for (unsigned i=0;i<Dumps.size();i++) {
        Dumps[i]->final_dump();
    }
    ///////////////////////////////////////////////
    // dump restart file for last configuration
    DumpLast->final_dump();
    ///////////////////////////////////////////////
    std::cout << "PolyMC finished!" << std::endl;
    return true;
}

bool PolyMC::run(long long steps ,const std::string & run_name,bool dump,bool print) {
    return run(steps,MCSteps,run_name,dump,print);
}

bool PolyMC::run(long long steps,std::vector<MCStep*> run_MCSteps ,const std::string & run_name,bool dump,bool print) {

    if (print) {
        std::cout << std::endl;
        std::cout << "######################################" << std::endl;
        std::cout << "######### Running Simulation #########" << std::endl;
        std::cout << "######################################" << std::endl;
        std::cout << std::endl;
    }

    start_timers();
    for (unsigned i=0;i<run_MCSteps.size();i++) {
        run_MCSteps[i]->reset_count();
    }
    if (check_link && EV_active) {
        set_link_backup();
    }
    if (EV_active) {
        EV->set_current_as_backup();
    }
    if (ES_active) {
        ES->set_current_as_backup(true);
    }
    chain->recal_energy();

    ///////////////////////////////////////////////
    ////////////// MAIN RUN LOOP //////////////////
    ///////////////////////////////////////////////

    long long step = 0;
    while (step<steps) {
        step++;

        if (step%print_every==0 && print) {
            print_state(step,steps,run_name,run_MCSteps);
        }

        if (step%check_consistency_every==0) {
            check_chain_consistency(step);
            if (check_link && EV_active) {
                check_link_consistency(step);
            }
        }
        for (unsigned i=0;i<run_MCSteps.size();i++) {
            run_MCSteps[i]->MC();
        }
        if (dump) {
            for (unsigned i=0;i<Dumps.size();i++) {
                Dumps[i]->dump();
            }
        }
    }
    check_chain_consistency(step);
    if (check_link && EV_active) {
        check_link_consistency(step);
    }
    return true;
}

void PolyMC::init_external_run() {
    external_run_steps = 0;
    start_timers();
    for (unsigned i=0;i<MCSteps.size();i++) {
        MCSteps[i]->reset_count();
    }
    if (check_link && EV_active) {
        set_link_backup();
    }
    if (EV_active) {
        EV->set_current_as_backup();
    }
    if (ES_active) {
        ES->set_current_as_backup(true);
    }
    chain->recal_energy();
}

bool PolyMC::external_run(long long steps ,const std::string & run_name,bool dump,bool print,bool recal_energy, bool set_backups, bool reset_count) {

    if (print) {
        std::cout << std::endl;
        std::cout << "######################################" << std::endl;
        std::cout << "######### Running Simulation #########" << std::endl;
        std::cout << "######################################" << std::endl;
        std::cout << std::endl;
    }

    if (reset_count) {
        external_run_steps = 0;
        start_timers();
        for (unsigned i=0;i<MCSteps.size();i++) {
            MCSteps[i]->reset_count();
        }
    }
    if (set_backups) {
        set_all_backups();
    }
    if (recal_energy) {
        chain->recal_energy();
    }

    ///////////////////////////////////////////////
    ////////////// MAIN RUN LOOP //////////////////
    ///////////////////////////////////////////////

    long long step = 0;
    while (step<steps) {
        step++;
        external_run_steps++;
//        if (external_run_steps%print_every==0 && print) {
//            print_state(step,steps,run_name,MCSteps);
//        }

        if (external_run_steps%check_consistency_every==0) {
//            std::cout << "checking chain consistency" << std::endl;
            check_chain_consistency(external_run_steps);
            if (check_link && EV_active) {
                check_link_consistency(external_run_steps);
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
//    check_chain_consistency(step);
//    if (check_link && EV_active) {
//        check_link_consistency(step);
//    }
    return true;
}



bool   PolyMC::exchange_configs(PolyMC * other_polymc) {
    arma::mat   temp_bp_pos = *chain->get_bp_pos();
    arma::cube  temp_triads = *chain->get_triads();
    double      temp_dLK    = chain->get_dLK();

    Chain * och = other_polymc->chain;

    chain->set_config(och->get_bp_pos(),och->get_triads(),och->topology_closed());
    chain->set_dLK(och->get_dLK());
    if (EV_active) {
        EV->set_backup_conf(och->get_bp_pos(),och->get_triads());
    }
    if (ES_active) {
        ES->set_current_as_backup(true);
    }
    chain->recal_energy();

    och->set_config(&temp_bp_pos,&temp_triads,chain->topology_closed());
    och->set_dLK(temp_dLK);
    if (EV_active) {
        other_polymc->EV->set_backup_conf(och->get_bp_pos(),och->get_triads());
    }
    if (ES_active) {
        other_polymc->ES->set_current_as_backup(true);
    }
    och->recal_energy();

    if (check_link && EV_active) {
        set_link_backup();
        other_polymc->set_link_backup();
    }
    return true;
}

double PolyMC::get_energy() {
    return chain->extract_energy();
}

Chain * PolyMC::get_chain() {
    return chain;
}


////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////

void PolyMC::init_general() {

    num_bp      = InputChoice_get_single<unsigned>    ("num_bp",input,argv,num_bp);
    num_bp      = InputChoice_get_single<unsigned>    ("nbp",input,argv,num_bp);
    steps       = InputChoice_get_single<long long>   ("steps",input,argv,steps);
    equi        = InputChoice_get_single<long long>   ("equi",input,argv,equi);
    temp        = InputChoice_get_single<double>      ("T",input,argv,temp);
    temp        = InputChoice_get_single<double>      ("temp",input,argv,temp);
    torque      = InputChoice_get_single<double>      ("torque",input,argv,force);

    std::vector<std::string> input_force_keys = {"force","f"};
    force       = InputChoice_get_single<double>      (input_force_keys,input,argv,force);

    closure_force  = InputChoice_get_single<double>         ("closure_force",input,argv,closure_force);
    closure_angularstiff = InputChoice_get_single<double>   ("closure_angularstiff",input,argv,closure_angularstiff);

    std::vector<std::string> input_sigma_keys = {"sig","sigma"};
    sigma       = InputChoice_get_single<double>      (input_sigma_keys,input,argv,sigma);

    std::vector<std::string> input_Lk0_from_static_keys = {"Lk0_from_static","Lk0_static"};
    Lk0_from_static = InputChoice_get_single<bool>        (input_Lk0_from_static_keys,input,argv,Lk0_from_static);

    std::vector<std::string> input_hel_rep_len_keys = {"hel_rep_len","hel_rep","helical_repeat","helical_repeat_length"};
    hel_rep_len = InputChoice_get_single<double>      (input_hel_rep_len_keys,input,argv,hel_rep_len);
    stateswitch_size      = InputChoice_get_single<unsigned>    ("stateswitch_size",input,argv,stateswitch_size);
    stateswitch_size      = InputChoice_get_single<unsigned>    ("ss_size",input,argv,stateswitch_size);
    stateswitch_num_moves      = InputChoice_get_single<unsigned>    ("stateswitch_num_moves",input,argv,stateswitch_num_moves);
    stateswitch_num_moves      = InputChoice_get_single<unsigned>    ("ss_num",input,argv,stateswitch_num_moves);

    geninfile.add_entry(GENINFILE_PARAMS,"steps",steps);
    geninfile.add_entry(GENINFILE_PARAMS,"equi",equi);
    geninfile.add_entry(GENINFILE_PARAMS,"num_bp",num_bp);
    geninfile.add_entry(GENINFILE_PARAMS,"T",temp);
    geninfile.add_entry(GENINFILE_PARAMS,"sigma",sigma);
    geninfile.add_entry(GENINFILE_PARAMS,"force",force);
    geninfile.add_entry(GENINFILE_PARAMS,"torque",torque);
    geninfile.add_entry(GENINFILE_PARAMS,"helical_repeat",hel_rep_len);
    geninfile.add_entry(GENINFILE_PARAMS,"Lk0_from_static",Lk0_from_static);
    geninfile.add_entry(GENINFILE_PARAMS,"ss_size",stateswitch_size);
    geninfile.add_entry(GENINFILE_PARAMS,"ss_num",stateswitch_num_moves);

    geninfile.add_entry(GENINFILE_PARAMS,"closure_force",closure_force);
    geninfile.add_entry(GENINFILE_PARAMS,"closure_angularstiff",closure_angularstiff);

    std::vector<std::string> input_IDB_keys = {"IDB_fn","IDB","idb"};
    IDB_fn      = InputChoice_get_single<std::string> (input_IDB_keys,input,argv,IDB_fn);
    subtract_T0 = InputChoice_get_single<bool>        ("subtract_T0",input,argv,subtract_T0);
    seq_fn = "";
    seq_fn = InputChoice_get_single<std::string> ("seqfn",input,argv,seq_fn);
    seq_fn = InputChoice_get_single<std::string> ("seq_fn",input,argv,seq_fn);
    seq_fn = InputChoice_get_single<std::string> ("sequence",input,argv,seq_fn);
//    seq_fn = input->get_single_val("seq_fn",seq_fn);
//    seq_fn = parse_arg(seq_fn, "-seq_fn",argv);

    std::vector<std::string> input_EV_keys = {"EV","EV_rad","ExVol"};
    EV_rad              = InputChoice_get_single<double>    (input_EV_keys,input,argv,EV_rad);
    EV_check_crossings  = InputChoice_get_single<bool>      ("EV_check_crossings",input,argv,EV_check_crossings);
    EV_force_repulsion_plane = InputChoice_get_single<bool> ("EV_force_repulsion_plane",input,argv,EV_force_repulsion_plane);

    geninfile.add_entry(GENINFILE_INTERACTIONS,"IDB",IDB_fn);
    geninfile.add_entry(GENINFILE_INTERACTIONS,"sequence",seq_fn);
    geninfile.add_entry(GENINFILE_INTERACTIONS,"subtract_T0",subtract_T0);
    geninfile.add_entry(GENINFILE_INTERACTIONS,"EV",EV_rad);
    geninfile.add_entry(GENINFILE_INTERACTIONS,"EV_check_crossings",EV_check_crossings);


    print_every       = InputChoice_get_single<int>   ("print_every",input,argv,print_every);
    print_link_info   = InputChoice_get_single<bool>  ("print_link_info",input,argv,print_link_info);
    print_elstat_info = InputChoice_get_single<bool>  ("print_elstat_info",input,argv,print_link_info);

    geninfile.add_entry(GENINFILE_OUTPUT,"print_every",print_every);
    geninfile.add_entry(GENINFILE_OUTPUT,"print_link_info",print_link_info);
    geninfile.add_entry(GENINFILE_OUTPUT,"print_elstat_info",print_elstat_info);


    use_cluster_twist = InputChoice_get_single<bool>  ("use_cluster_twist",input,argv,use_cluster_twist);
    num_twist         = InputChoice_get_single<unsigned> ("num_twist",input,argv,num_twist);
    check_link        = InputChoice_get_single<bool>  ("check_link",input,argv,check_link);
    check_consistency = InputChoice_get_single<bool>  ("check_consistency",input,argv,check_consistency);
    check_consistency_every = InputChoice_get_single<int>   ("check_consistency_every",input,argv,check_consistency_every);
    if (check_consistency_every==0) {
        check_consistency_every=1e15;
    }

    geninfile.add_entry(GENINFILE_SIMSETUP,"use_cluster_twist",use_cluster_twist);
    geninfile.add_entry(GENINFILE_SIMSETUP,"num_twist",num_twist);
    geninfile.add_entry(GENINFILE_SIMSETUP,"check_link",check_link);
    geninfile.add_entry(GENINFILE_SIMSETUP,"check_consistency",check_consistency);
    geninfile.add_entry(GENINFILE_SIMSETUP,"check_consistency_every",check_consistency_every);


    copy_input        = InputChoice_get_single<bool>  ("copy_input",input,argv,copy_input);
    geninfile.add_entry(GENINFILE_DUMPS_SETUP,"copy_input",copy_input);


    restart_file     = InputChoice_get_single<std::string>  ("restart",input,argv,restart_file);
    restart_snapshot = InputChoice_get_single<int>          ("restart_snapshot",input,argv,restart_snapshot);
    restart_set_link = InputChoice_get_single<bool>         ("restart_set_link",input,argv,restart_set_link);
    if (restart_file != "") {
        geninfile.add_entry(GENINFILE_SIMSETUP,"restart",restart_file);
        geninfile.add_entry(GENINFILE_SIMSETUP,"restart_snapshot",restart_snapshot);
        geninfile.add_entry(GENINFILE_SIMSETUP,"restart_set_link",restart_set_link);
    }



    constraint_fn    = InputChoice_get_single<std::string>  ("constraint_fn",input,argv,constraint_fn);
    if (constraint_fn != "") {
        geninfile.add_entry(GENINFILE_INTERACTIONS,"constraint_fn",constraint_fn);
    }

    GenForce_fn      = InputChoice_get_single<std::string>  ("GenForce_fn",input,argv,GenForce_fn);
    GenForce_fn      = InputChoice_get_single<std::string>  ("genForce_fn",input,argv,GenForce_fn);
    GenForce_fn      = InputChoice_get_single<std::string>  ("genforce_fn",input,argv,GenForce_fn);
    if (GenForce_fn != "") {
        geninfile.add_entry(GENINFILE_INTERACTIONS,"GenForce_fn",GenForce_fn);
    }

    ElStat_fn       = InputChoice_get_single<std::string>  ("ElStat_fn",input,argv,ElStat_fn);
    ElStat_fn       = InputChoice_get_single<std::string>  ("elStat_fn",input,argv,ElStat_fn);
    ElStat_fn       = InputChoice_get_single<std::string>  ("elstat_fn",input,argv,ElStat_fn);
    ElStat_fn       = InputChoice_get_single<std::string>  ("Electrostatics_fn",input,argv,ElStat_fn);
    if (ElStat_fn != "") {
        geninfile.add_entry(GENINFILE_INTERACTIONS,"ElStat_fn",ElStat_fn);
    }

    pair_interactions_fn      = InputChoice_get_single<std::string>  ("set_pairstyles",input,argv,pair_interactions_fn);
    initialize_atom_types_fn  = InputChoice_get_single<std::string>  ("set_atoms",input,argv,initialize_atom_types_fn);
    if (pair_interactions_fn != "") {
        geninfile.add_entry(GENINFILE_INTERACTIONS,"set_pairstyles",pair_interactions_fn);
    }
    if (initialize_atom_types_fn != "") {
        geninfile.add_entry(GENINFILE_INTERACTIONS,"set_atoms",initialize_atom_types_fn);
    }

    std::cout << " num_bp           = " << num_bp << std::endl;
    std::cout << " equi             = " << equi << std::endl;
    std::cout << " steps            = " << steps << std::endl;
    std::cout << " IDB              = " << IDB_fn << std::endl;
    std::cout << " temp             = " << temp << std::endl;
    std::cout << " torque           = " << torque << std::endl;
    std::cout << " force            = " << force << std::endl;
    std::cout << " sigma            = " << sigma << std::endl;
    std::cout << " subtract_T0      = " << subtract_T0 << std::endl;
    std::cout << " Lk0_from_static  = " << Lk0_from_static << std::endl;
    std::cout << " print_every      = " << print_every << std::endl;

    if (closure_force != 0 || closure_angularstiff != 0) {
        std::cout << " closure_force        = " << print_every << std::endl;
        std::cout << " closure_angularstiff = " << print_every << std::endl;
    }

    if (check_consistency_every) {
        std::cout << std::endl << " checking consistency every " << check_consistency_every << " steps." << std::endl;
    }

}

void PolyMC::init_seed() {
    seedstr = InputChoice_get_single<std::string>  ("seed",input,argv,seedstr);
    seedseq = seedstr2seedseq(seedstr);
    seedstr = seedseq2seedstr(seedseq);

//    if (seedstr == "-1") {
//        std::random_device  rd{};
//        for (int i=0;i<SEED_SEQ_LEN;i++) {
//            seedseq.push_back(rd());
//        }
//        std::ostringstream ss;
//        ss << seedseq[0];
//        seedstr = ss.str();
//        for (int i=1;i<SEED_SEQ_LEN;i++) {
//            std::ostringstream ss1;
//            ss1 << seedseq[i];
//            seedstr += "-"+ss1.str();
//        }
//    }
//    else {
//        size_t pos = 0;
//        std::string num;
//        std::string cp_seedstr = seedstr;
//        while ((pos = cp_seedstr.find("-")) != std::string::npos) {
//            num = cp_seedstr.substr(0, pos);
//            seedseq.push_back(std::stoll(num));
//            cp_seedstr.erase(0, pos + 1);
//        }
//        seedseq.push_back(std::stoll(cp_seedstr));
//    }
    geninfile.add_entry(GENINFILE_SEEDS,"seed",seedstr);
}

void PolyMC::init_geninfile() {
    geninfile.add_category(GENINFILE_SIMSETUP);
    geninfile.add_category(GENINFILE_SEEDS);
    geninfile.add_category(GENINFILE_INTERACTIONS);
    geninfile.add_category(GENINFILE_PARAMS);
    geninfile.add_category(GENINFILE_OUTPUT);
    geninfile.add_category(GENINFILE_DUMPS_SETUP);
    geninfile.add_category(GENINFILE_DUMPS);
}

void PolyMC::init_dump() {

    if (input_given) {
        dump_dir = input->get_single_val<std::string>("dump_dir",dump_dir);
    }
    dump_dir = parse_arg(dump_dir, "-dir"      ,argv);
    dump_dir = parse_arg(dump_dir, "-dump_dir" ,argv);
    geninfile.add_entry(GENINFILE_DUMPS_SETUP,"dump_dir",dump_dir);

    Dumps = init_dump_cmdargs(argv,&geninfile,chain,mode,EV_rad,inputfn);
    /*
        Dump Energy
    */
    /*
        For now I placed this one here to avoid having to pass the Elstat pointer as an argument
    */

    int num_bp  = chain->get_num_bp();
    std::string dump_dir = "";

    bool inputfile_given=false;
    InputRead input(inputfn);
    if (input.filefound()){
        dump_dir  = input.get_single_val<std::string>("dump_dir",dump_dir);
        inputfile_given=true;
    }
    dump_dir    = parse_arg(dump_dir, "-dir"     ,argv);
    dump_dir    = parse_arg(dump_dir, "-dump_dir",argv);

    bool append_dumps = false;
    append_dumps = parse_flag(append_dumps,"-app", argv);
    append_dumps = InputChoice_get_single<bool>        ("app"     ,&input,argv,append_dumps);
    append_dumps = InputChoice_get_single<bool>        ("append"  ,&input,argv,append_dumps);
    append_dumps = InputChoice_get_single<bool>        ("append"  ,&input,argv,append_dumps);

    geninfile.add_entry(GENINFILE_DUMPS_SETUP,"app",append_dumps);
    geninfile.add_newline(GENINFILE_DUMPS_SETUP);

    if (ES_active) {

        int    ESEnergy_dump_every      = 0;
        std::string ESEnergy_filename   = dump_dir+".esen";
        bool   ESENERGY_kt_units        = false;
        bool   ESENERGY_recal           = true;

        ESEnergy_dump_every  = InputChoice_get_single<int>         ("ESEn" ,&input,argv,ESEnergy_dump_every);
        ESEnergy_dump_every  = InputChoice_get_single<int>         ("esen" ,&input,argv,ESEnergy_dump_every);
        ESEnergy_filename    = InputChoice_get_single<std::string> ("ESEfn",&input,argv,ESEnergy_filename);
        ESEnergy_filename    = InputChoice_get_single<std::string> ("esefn",&input,argv,ESEnergy_filename);

        ESENERGY_kt_units    = InputChoice_get_single<bool>        ("ESEkt" ,&input,argv,ESENERGY_kt_units);
        ESENERGY_kt_units    = InputChoice_get_single<bool>        ("esekt" ,&input,argv,ESENERGY_kt_units);

        ESENERGY_recal       = InputChoice_get_single<bool>        ("ESErecal" ,&input,argv,ESENERGY_recal);
        ESENERGY_recal       = InputChoice_get_single<bool>        ("eserecal" ,&input,argv,ESENERGY_recal);

        if (ESEnergy_dump_every>0) {

            Dump_ESEnergy* DESEnergy = new Dump_ESEnergy(chain, ESEnergy_dump_every, ESEnergy_filename ,ES,ESENERGY_kt_units, ESENERGY_recal,append_dumps);
            Dumps.push_back(DESEnergy);

            geninfile.add_entry(GENINFILE_DUMPS,"ESEn" ,ESEnergy_dump_every);
    //        geninfile->add_entry(GENINFILE_DUMPS,"Efn",Energy_filename);
            geninfile.add_newline(GENINFILE_DUMPS);
        }
    }


    std::cout << std::endl << Dumps.size() << " Dumps initialized!" << std::endl << std::endl;

    Dump_Restart * nlast = new Dump_Restart(chain, 0, dump_dir+".last",false,true);
    DumpLast = nlast;
}

void PolyMC::init_ExVol() {

//    if (EV_rad <= 0) {
//        EV_active=false;
//        return;
//    }

    if (EV_rad >= chain->get_disc_len() || (EV_force_repulsion_plane && EV_rad == 0)) {
        EV_active=true;
        ExVol * exvol = new ExVol(chain,EV_rad,EV_check_crossings);
        EV = exvol;

        if (EV_repulsion_plane) {

            if (tweezer_boundary_front != "surface" && tweezer_boundary_front != "line" && tweezer_boundary_front != "none") {
                std::cout << "Error: Unknown tweezer_boundary_front setting '" << tweezer_boundary_front << "'!" << std::endl;
                std::cout << " Valid settings: " << std::endl;
                std::cout << "  - 'surface'" << std::endl;
                std::cout << "       -> Repulsion planes" << std::endl;
                std::cout << "  - 'line'" << std::endl;
                std::cout << "       -> Ghost closure to infinity" << std::endl;
                std::cout << "  - 'none'" << std::endl;
                std::cout << "       -> No boundary closure" << std::endl;
                std::exit(0);
            }
            if (tweezer_boundary_back != "surface" && tweezer_boundary_back != "line" && tweezer_boundary_back != "none") {
                std::cout << "Error: Unknown tweezer_boundary_back setting '" << tweezer_boundary_back << "'!" << std::endl;
                std::cout << " Valid settings: " << std::endl;
                std::cout << "  - 'surface'" << std::endl;
                std::cout << "       -> Repulsion planes" << std::endl;
                std::cout << "  - 'line'" << std::endl;
                std::cout << "       -> Ghost closure to infinity" << std::endl;
                std::cout << "  - 'none'" << std::endl;
                std::cout << "       -> No boundary closure" << std::endl;
                std::exit(0);
            }

            bf_surface  = false;
            bb_surface  = false;
            bf_line     = false;
            bb_line     = false;

            if (tweezer_boundary_front == "surface") {
                bf_surface = true;
            }
            if (tweezer_boundary_back  == "surface") {
                bb_surface = true;
            }
            if (tweezer_boundary_front == "line") {
                bf_line = true;
            }
            if (tweezer_boundary_back  == "line") {
                bb_line = true;
            }

            if (bf_surface || bb_surface) {
                EV_repulsion_plane = true;
                EV->set_repulsion_plane(bf_surface,bb_surface);
            }
            else {
                EV_repulsion_plane = false;
            }

            if (bf_line || bb_line) {
                EV_line_closure = true;
                EV->set_line_closure(bf_line,bb_line);
            }
            else {
                EV_line_closure = false;
            }
        }

//        EV->set_repulsion_plane(EV_repulsion_plane);
//        EV->set_line_closure(EV_line_closure);
//        EV->set_repulsion_plane(false,false);
//        EV->set_line_closure(true,true);


        for (unsigned i=0;i<MCSteps.size();i++) {
            MCSteps[i]->set_excluded_volume(EV);
        }
    }
    else {
        if (EV_rad > 0) {

            std::cout << "Excluded Volume radius (" << EV_rad << ") is smaller than discretization length. Simulation terminated!" << std::endl;
            std::exit(0);
        }
        else {
            EV_active=false;
        }
    }
    std::cout << std::endl;
}

void PolyMC::init_ElStat(const std::string & ElStat_fn) {
    /*
        TODO: Check for double entries
    */

    if (ElStat_fn == "") {
        print_elstat_info = false;
        return;
    }


    /*
        Read Electrostatics
    */
    InputRead * input = new InputRead(ElStat_fn);
    std::string type            = input->get_single_val<std::string> ("type");
    std::vector<double> params  = input->get_single_vec("params");
    double integral_dx          = input->get_single_val<double> ("integral_dx");

    double neighbor_skip_dist   = input->get_single_val<double> ("neighbor_skip_dist");
    double cutoff_dist          = input->get_single_val<double> ("cutoff_dist");

    bool tabulate               = input->get_single_val<bool>   ("tabulate");
    bool use_costheta           = input->get_single_val<bool>   ("use_costheta");
    int  table_elems            = input->get_single_val<double> ("table_elements");
    table_elems                 = input->get_single_val<double> ("table_elems",table_elems);
    double memory_limit         = input->get_single_val<double> ("memory_limit",1000);

    double rho_max,rho_min;
    rho_max = cutoff_dist + chain->get_disc_len();
    ES_active = init_electrostatic_potential(type,params,integral_dx,rho_max);
    if (EV_active) {
        rho_min = EV_rad;
    }
    else {
        rho_min = 0;
    }

    if (tabulate) {
        ES_Potential->tabulate_interactions(table_elems,rho_min,rho_max,memory_limit,use_costheta);
    }
    ElStat * elstat = new ElStat(chain,ES_Potential,rho_max,rho_min,neighbor_skip_dist );
    ES = elstat;
    ES->set_backup_pos(EV->get_bp_pos_backup(),EV->get_triads_backup());

    for (unsigned i=0;i<MCSteps.size();i++) {
        MCSteps[i]->set_electrostatics(ES);
    }
}





void PolyMC::init_constraints(const std::string & constraint_fn, bool append) {
    /*
        TODO: check for duplicates and remove them
    */

    if (constraint_fn == "") return;

    if (!append) {
        for (unsigned i=0;i<constraints.size();i++) {
            delete constraints[i];
        }
        constraints.clear();
    }

    InputRead * constraint_read = new InputRead(constraint_fn);

    /*
        Init Position Constraints
    */
    unsigned fix_pos_num = constraint_read->get_single_num_instances("fix_pos");
    for (int i=0;i<fix_pos_num;i++) {
        std::vector<int> fix_group = constraint_read->get_single_intvec(i,"fix_pos");
        if (EV_active) {
            Constr_Position * fix_pos = new Constr_Position(chain,fix_group, EV->get_bp_pos_backup(), EV->get_triads_backup());
            constraints.push_back(fix_pos);
        }
        else {
            Constr_Position * fix_pos = new Constr_Position(chain,fix_group);
            constraints.push_back(fix_pos);
        }
    }

    /*
        Init Orientation Constraints
    */
    unsigned fix_ori_num = constraint_read->get_single_num_instances("fix_ori");
    for (int i=0;i<fix_ori_num;i++) {
        std::vector<int> fix_group = constraint_read->get_single_intvec(i,"fix_ori");
        if (EV_active) {
            Constr_Orientation * fix_ori = new Constr_Orientation(chain,fix_group, EV->get_bp_pos_backup(), EV->get_triads_backup());
            constraints.push_back(fix_ori);
        }
        else {
            Constr_Orientation * fix_ori = new Constr_Orientation(chain,fix_group);
            constraints.push_back(fix_ori);
        }
    }

    if (constraints.size() > 0) {
        constraints_active = true;
        for (unsigned i=0;i<MCSteps.size();i++) {
            MCSteps[i]->set_constraints(constraints);
        }
    }
    std::cout << constraints.size() << " Constraints found!" << std::endl;
}


void PolyMC::init_constraints() {

    if (input->contains_multi("Constraints:")) {

        MultiLineInstance * ML = input->get_multiline("Constraints:")->get_instance(0);
        std::vector<SingleLineInstance> * singlelines = ML->get_all_singlelines();

        for (unsigned i=0;i<singlelines->size();i++) {

            if ((*singlelines)[i].identifier == "fix_position" || (*singlelines)[i].identifier == "fix_pos") {
                std::vector<int> fix_group = (*singlelines)[i].get_intvec();
                if (EV_active) {
                    Constr_Position * fix_pos = new Constr_Position(chain,fix_group, EV->get_bp_pos_backup(), EV->get_triads_backup());
                    constraints.push_back(fix_pos);
                }
                else {
                    Constr_Position * fix_pos = new Constr_Position(chain,fix_group);
                    constraints.push_back(fix_pos);
                }
            }

            if ((*singlelines)[i].identifier == "fix_orientation" || (*singlelines)[i].identifier == "fix_ori") {
                std::vector<int> fix_group = (*singlelines)[i].get_intvec();
                if (EV_active) {
                    Constr_Orientation * fix_ori = new Constr_Orientation(chain,fix_group, EV->get_bp_pos_backup(), EV->get_triads_backup());
                    constraints.push_back(fix_ori);
                }
                else {
                    Constr_Orientation * fix_ori = new Constr_Orientation(chain,fix_group);
                    constraints.push_back(fix_ori);
                }
            }

        }
        if (constraints.size() > 0) {
            constraints_active = true;
            for (unsigned i=0;i<MCSteps.size();i++) {
                MCSteps[i]->set_constraints(constraints);
            }
//            std::exit(1);
        }
    }
    std::cout << constraints.size() << " Constraints found!" << std::endl;
}

void PolyMC::init_GenForce(const std::string & GenForce_fn) {

    /*
        TODO: Check for double entries
    */

    if (GenForce_fn == "") return;

    InputRead * GenForce_read = new InputRead(GenForce_fn);
    /*
        Read Generalized Force
    */
    std::vector<std::vector<double>> list_GenForce;
    unsigned num = GenForce_read->get_single_num_instances("GenF");
    for (int i=0;i<num;i++) {
        std::vector<double> gf = GenForce_read->get_single_vec(i,"GenF");
        if (gf.size() >= 4) {
            list_GenForce.push_back(gf);
        }
        else {
            std::cout << "Warning: Insufficient entries provided for GenForce! 4 entries are required, " << gf.size() << " given." << std::endl;
        }
    }
    for (unsigned i=0;i<list_GenForce.size();i++) {
        int bps_id = list_GenForce[i][0];
        arma::colvec gf_vec = {list_GenForce[i][1],list_GenForce[i][2],list_GenForce[i][3]};

        chain->get_BPS()->at(bps_id)->set_generalized_force(gf_vec);
    }
}



void PolyMC::init_pair_interactions() {

    if (pair_interactions_fn!="" && initialize_atom_types_fn!="") {
        pair_interactions_active = true;
        unbound = new Unbound(chain);
        unbound->read_inputfile(pair_interactions_fn);
        unbound->init_monomer_type_assignments_from_file(initialize_atom_types_fn);
        unbound->set_backup();

        for (unsigned i=0;i<MCSteps.size();i++) {
            MCSteps[i]->set_unbound(unbound);
        }
    }
}

void PolyMC::set_all_backups() {
    if (EV_active) {
        EV->set_current_as_backup();
    }
    if (ES_active) {
        ES->set_current_as_backup(true);
    }
    if (check_link && EV_active) {
        set_link_backup();
    }
}


void PolyMC::init_stateswitch_moves() {
    if (stateswitch_size == 0 || stateswitch_num_moves == 0) {
        return;
    }
    for (unsigned i=0;i<chain->get_BPS()->size();i++) {
        if (!chain->get_BPS()->at(i)->is_bubble()) {
            stateswitch_size = 0;
            return;
        }
    }
    if (stateswitch_size > 0) {
        std::cout << " Adding Stateswitch to MC Moves .. " << std::endl;
        std::cout << stateswitch_size << std::endl;
        MCS_Stateswitch*  sw  = new MCS_Stateswitch(chain,seedseq,0,larger(0,stateswitch_size-1));
        for (unsigned jj=0;jj<stateswitch_num_moves;jj++) {
            MCSteps.push_back(sw);
        }
    }
}


bool PolyMC::check_chain_consistency(long long step) {
    bool consistent = chain->check_energy_consistency();
    if (!consistent) {
        std::cout << "Energy inconsistent!" << std::endl;
        std::string energyviolatedfn = dump_dir + ".energyinconsistent";
        std::ofstream ofstr;
        ofstr.open(energyviolatedfn, std::ofstream::out | std::ofstream::app);
        ofstr << step << std::endl;
        ofstr.close();
        // #ifdef POLYMC_TERMINATE_INCONSISTENT_ENERGY
        // std::exit(0);
        // #endif
        chain->recal_energy();
    }
    if (!chain->config_consistent()) {
        std::cout << "chain positions inconsistent!" << std::endl;

        std::string posviolatedfn = dump_dir + ".posinconsistent";
        std::ofstream ofstr;
        ofstr.open(posviolatedfn, std::ofstream::out | std::ofstream::app);
        ofstr << step << std::endl;
        ofstr.close();
        // #ifdef POLYMC_TERMINATE_INCONSISTENT_CONFIG
        // std::exit(0);
        // #endif
        chain->restore_consistency();
        if (EV_active) {
            if (EV->check_overlap()) {
                std::cout << "Consistency check lead to overlap!" << std::endl;
                EV->revert_to_backup();
            }
            else {
                EV->set_current_as_backup();
            }
        }
        if (ES_active) {
            ES->set_current_as_backup(true);
        }
        if (check_link && EV_active) {
            set_link_backup();
        }
        chain->recal_energy();
    }
    return true;
}


void PolyMC::set_link_backup() {
    check_link_backup_triads = *chain->get_triads();
    check_link_backup_bp_pos = *chain->get_bp_pos();
    check_link_backup_dLK    =  chain->get_dLK();
}


bool PolyMC::check_link_consistency(long long step) {

    double Wr,Tw,dLK;
//    Wr  = chain->cal_langowski_writhe_1a();

    arma::mat writhe_elem;
    chain->langowski_writhe_elements(&writhe_elem, false);
    Wr = arma::accu(writhe_elem);

    Tw  = chain->cal_twist(0,num_bp);
    dLK = chain->get_dLK();

    if (std::abs(dLK-Tw-Wr)> LINK_CONSERVATION_VIOLATION_THRESHOLD) {
        std::cout << "  LK = " << Tw+Wr << " ( " << chain->get_dLK() << ")" << std::endl;
        std::cout << "Linking Number violated... reverting to backup configuration." << std::endl;
        chain->set_config(&check_link_backup_bp_pos, &check_link_backup_triads, false, false);
        chain->set_dLK(check_link_backup_dLK);
        EV->set_current_as_backup();
        if (ES_active) {
            ES->set_current_as_backup(true);
        }

        std::string lkviolatedfn = dump_dir + ".linkviolated";
        std::ofstream ofstr;
        ofstr.open(lkviolatedfn, std::ofstream::out | std::ofstream::app);
        ofstr << step << " " << dLK << " " << Tw+Wr << " " << Wr << " " << Tw << std::endl;
        ofstr.close();
        return false;
    }
    else {
        check_link_backup_triads = *chain->get_triads();
        check_link_backup_bp_pos = *chain->get_bp_pos();
        check_link_backup_dLK    =  chain->get_dLK();
        return true;
    }
}

double PolyMC::cal_dlk() {
    double Wr,Tw;
//    Wr  = chain->cal_langowski_writhe_1a();
    arma::mat writhe_elem;
    chain->langowski_writhe_elements(&writhe_elem, false);
    Wr = arma::accu(writhe_elem);
    Tw  = chain->cal_twist(0,num_bp);
    // std::cout << "  LK = " << Tw+Wr << " ( " << chain->get_dLK() << ")" << std::endl;
    return Tw+Wr;
}



void PolyMC::start_timers() {
    timer_start  = std::chrono::high_resolution_clock::now();
    timer_finish = std::chrono::high_resolution_clock::now();
    timer_ref    = std::chrono::high_resolution_clock::now();
}

void PolyMC::print_state(long long step,long long steps,const std::string & run_name,const std::vector<MCStep*> & run_MCSteps) {

    double Wr,Tw;
    if (print_link_info) {

        if (EV_line_closure) {
            int frontadd = 0;
            int backadd  = 0;

            if (bf_line) {
                frontadd = int(num_bp/2);
            }
            if (bb_line) {
                backadd = int(num_bp/2);
            }
            Wr = chain->Chain::cal_langowski_writhe_1a(frontadd, backadd);
        }
        else {
            Wr = chain->cal_langowski_writhe_1a();
        }
        Tw = chain->cal_twist(0,num_bp);
    }

    std::cout << std::endl;
    std::cout << "################################################" << std::endl;
    std::cout << " " << run_name << ": --- step " << step << " / " << steps << " --- " << std::endl << std::endl;
    std::cout << " MCSteps:" << std::endl;
    print_acceptance_rates(run_MCSteps);


    if (mode==POLYMC_MODE_TWEEZER && tweezer_setup_id>=0) {
        std::cout << std::endl;
        std::cout << " Tweezer Info:" << std::endl;
        std::cout << "  Current dLK = " << chain->get_dLK() << std::endl;
        if (tweezer_setup_id==TWEEZER_SETUP_ID_TORSIONAL_TRAP) {
            std::cout << "  Mean Torque = " << chain->get_torque_measure_accu(false)/(2*M_PI) << std::endl;
        }
    }
    if (print_elstat_info) {
        std::cout << std::endl;
        std::cout << " Electrostatics Info:" << std::endl;
        double Ees = ES->get_current_energy(true);
        std::cout << "  E_es = " << Ees << " (kT)" << std::endl;

    }
    if (print_link_info) {
        std::cout << std::endl;
        std::cout << " Link Info:" << std::endl;
        std::cout << "  Wr = " << Wr << std::endl;
        std::cout << "  Tw = " << Tw << std::endl;
        std::cout << "  LK = " << Tw+Wr << " ( " << chain->get_dLK() << ")" << std::endl;
    }
    std::cout << std::endl;
    std::cout << " Time Info:" << std::endl;
    timer_finish   = std::chrono::high_resolution_clock::now();
    timer_elapsed  = timer_finish - timer_ref;
    std::cout << "  elapsed time:             " << timer_elapsed.count() << "s\n";
    timer_ref = std::chrono::high_resolution_clock::now();
    timer_elapsed = timer_finish-timer_start;
    std::cout << "  estimated to finished in: " << time_remaining((timer_elapsed.count())/step*(steps-step)) << "\n";

    // std::cout << (*chain->get_states()).t() << std::endl;
    // std::cout << arma::sum((*chain->get_states()))- chain->get_states()->size() << std::endl;

////    // TESTING !!!!
//    int dim = 1;
//    std::cout << " ---------" << std::endl;
//    std::cout << " ---------" << std::endl;
//    std::cout << chain->get_triads()->slice(19).col(dim) << std::endl;
//    std::cout << chain->get_triads()->slice(20).col(dim) << std::endl;
//    std::cout << chain->get_triads()->slice(21).col(dim) << std::endl;
//
//    int dim = 1;
//    std::cout << " ---------" << std::endl;
//    std::cout << " ---------" << std::endl;
//    std::cout << chain->get_triads()->slice(2).col(dim) << std::endl;
//    std::cout << chain->get_triads()->slice(3).col(dim) << std::endl;
//    std::cout << chain->get_triads()->slice(4).col(dim) << std::endl;
//
//    std::cout << " ---------" << std::endl;
//    std::cout << chain->get_triads()->slice(400).col(dim) << std::endl;
//    std::cout << chain->get_triads()->slice(401).col(dim) << std::endl;
//    std::cout << chain->get_triads()->slice(402).col(dim) << std::endl;


//    std::cout << chain->get_bp_pos()->col(4).col(dim) << std::endl;
//    std::cout << chain->get_bp_pos()->col(5).col(dim) << std::endl;
//    std::cout << chain->get_bp_pos()->col(6).col(dim) << std::endl;
//
//    std::cout << num_bp << std::endl;
//    std::cout << chain->get_bp_pos()->n_cols << std::endl;
//    std::cout << chain->get_disc_len() << std::endl;
    // !!!!!!!!!!!!!!!!!!!!!!!!
}


void PolyMC::print_acceptance_rates(const std::vector<MCStep*> & run_MCSteps) {
    for (unsigned i=0;i<run_MCSteps.size();i++) {
        run_MCSteps[i]->print_acceptance_rate();
    }
}


std::string PolyMC::time_remaining(int seconds) {
    int days    = seconds/86400;
    seconds     = seconds%86400;
    int hours   = seconds/3600;
    seconds     = seconds%3600;
    int minutes = seconds/60;
    seconds     = seconds%60;
    std::string str;
    if (days>0) {
        str += std::to_string(days) + "d ";
    }
    if (hours>0) {
        str += std::to_string(hours) + "h ";
    }
    if (minutes>0) {
        str += std::to_string(minutes) + "min ";
    }
    str += std::to_string(seconds) + "s ";
    return str;
}


void PolyMC::copy_input_files() {
//    std::string copy_inputfn = dump_dir+".in";
//    std::ifstream  src_in(inputfn,      std::ios::binary);
//    std::ofstream  dst_in(copy_inputfn, std::ios::binary);
//    dst_in << src_in.rdbuf();

    std::string gen_inputfn = dump_dir+".in";
    geninfile.generate_input_file(gen_inputfn);

    std::string copy_IDBfn   = dump_dir+".idb";
    std::ifstream  src_idb(IDB_fn,     std::ios::binary);
    std::ofstream  dst_idb(copy_IDBfn, std::ios::binary);
    dst_idb << src_idb.rdbuf();

    if (constraint_fn != "") {
        std::string    copy_constr_fn   = dump_dir+".constr";
        std::ifstream  src_constr(constraint_fn,  std::ios::binary);
        std::ofstream  dst_constr(copy_constr_fn, std::ios::binary);
        dst_constr << src_constr.rdbuf();
    }

    if (GenForce_fn != "") {
        std::string    copy_GenForce_fn   = dump_dir+".genForce";
        std::ifstream  src_gf(GenForce_fn,  std::ios::binary);
        std::ofstream  dst_gf(copy_GenForce_fn, std::ios::binary);
        dst_gf << src_gf.rdbuf();
    }

    if (ElStat_fn != "") {
        std::string    copy_ElStat_fn   = dump_dir+".elstat";
        std::ifstream  src_gf(ElStat_fn,  std::ios::binary);
        std::ofstream  dst_gf(copy_ElStat_fn, std::ios::binary);
        dst_gf << src_gf.rdbuf();
    }

    if (seq_fn != "") {
        std::string    copy_seq_fn   = dump_dir+".seq";
        std::ifstream  src_seq(seq_fn,  std::ios::binary);
        std::ofstream  dst_seq(copy_seq_fn, std::ios::binary);
        dst_seq << src_seq.rdbuf();
    }
}

bool PolyMC::generic_slow_wind(double dLK_from, double dLK_to, double dLK_step, long long steps_per_dLK_step) {

    int sgn_of_change = sgn(dLK_to-dLK_from);
    if (sgn(dLK_step) != sgn_of_change) {
        dLK_step = -dLK_step;
    }

    double set_dLK = dLK_from;
    chain->set_Delta_Lk(set_dLK);
    EV->set_current_as_backup();
    if (ES_active) {
        ES->set_current_as_backup(true);
    }

    if (check_link) {
        set_link_backup();
    }

    while (set_dLK*sgn_of_change < dLK_to*sgn_of_change ) {
        std::cout << "dLK set to " << set_dLK << std::endl;
        run(steps_per_dLK_step,MCSteps,"winding chain - target: "+to_string_with_precision(dLK_to,2)+" turns\n",slow_wind_dump);
        set_dLK += dLK_step;

        if (set_dLK*sgn_of_change > dLK_to*sgn_of_change) {
            set_dLK = dLK_to;
        }
        chain->set_Delta_Lk(set_dLK);
        EV->set_current_as_backup();
        if (ES_active) {
            ES->set_current_as_backup(true);
        }
    }
    std::cout << "dLK set to " << set_dLK << std::endl;
    run(steps_per_dLK_step,MCSteps,"winding chain",slow_wind_dump);

    double Wr,Tw;
    Wr  = chain->cal_langowski_writhe_1a();
    Tw  = chain->cal_twist(0,num_bp);
    std::cout << "  LK = " << Tw+Wr << " ( " << chain->get_dLK() << ")" << std::endl;

    return true;
//    for (unsigned i=0;i<wind_MCSteps.size();i++) {
//        delete wind_MCSteps[i];
//    }
}


