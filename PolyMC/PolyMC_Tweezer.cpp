#include "../PolyMC.h"

bool PolyMC::init_tweezer() {

    std::cout << std::endl;
    std::cout << "######################################" << std::endl;
    std::cout << "###### Loading Tweezer Protocol ######" << std::endl;
    std::cout << "######################################" << std::endl;
    std::cout << std::endl;

    closed                  = false;
    terminus_fixed_tangents = true;
    terminus_fixed_triads   = false;
    EV_repulsion_plane      = true;

    tweezer_boundary_front  = InputChoice_get_single<std::string>  ("tweezer_boundary_front",input,argv,tweezer_boundary_front);
    tweezer_boundary_back   = InputChoice_get_single<std::string>  ("tweezer_boundary_back" ,input,argv,tweezer_boundary_back);

    geninfile.add_entry(GENINFILE_SIMSETUP,"tweezer_boundary_front",tweezer_boundary_front);
    geninfile.add_entry(GENINFILE_SIMSETUP,"tweezer_boundary_back",tweezer_boundary_back);


//    EV_line_closure         = false;
//
//    EV_line_closure         = InputChoice_get_single<bool>  ("tweezer_boundary_line",input,argv,EV_line_closure);
//    if (EV_line_closure) EV_repulsion_plane  = false;
//
//    EV_repulsion_plane      = InputChoice_get_single<bool>  ("tweezer_boundary_surface",input,argv,EV_repulsion_plane);
//    if (EV_repulsion_plane) EV_line_closure = false;



    /*
        Tweezer Settings
    */
    tweezer_setup_id = -1;
    tweezer_setup        = InputChoice_get_single<std::string> ("tweezer_setup",input,argv,tweezer_setup);
    geninfile.add_entry(GENINFILE_SIMSETUP,"tweezer_setup",tweezer_setup);
    if (tweezer_setup==TWEEZER_SETUP_FREE) {
        tweezer_setup_id = TWEEZER_SETUP_ID_FREE;
    }
    if (tweezer_setup==TWEEZER_SETUP_LINKFIX) {
        tweezer_setup_id = TWEEZER_SETUP_ID_LINKFIX;
        terminus_fixed_triads   = true;
    }
    if (tweezer_setup==TWEEZER_SETUP_TORSIONAL_TRAP) {
        tweezer_setup_id = TWEEZER_SETUP_ID_TORSIONAL_TRAP;
    }
    if (tweezer_setup==TWEEZER_SETUP_TORQUE) {
        tweezer_setup_id = TWEEZER_SETUP_ID_TORQUE;
        sigma = 0;
        std::cout << "Warning: Tweezer setup Torque requires the initial supercoiling" << std::endl;
        std::cout << "density to be set to zero. (sigma was set to 0)" << std::endl;
    }

    /*
        For vids
    */
    bool dynamics = false;
    if (tweezer_setup==TWEEZER_SETUP_LINKFIX_DYNAMICS) {
        tweezer_setup_id        = TWEEZER_SETUP_ID_LINKFIX;
        terminus_fixed_triads   = true;
        dynamics                = true;
    }

    if (tweezer_setup_id<0) {
        std::cout << "Error: '" << tweezer_setup << "' is not a valid tweezer setup!" << std::endl;
        std::cout << " Valid setups:" << std::endl;
        std::cout << "   - " << TWEEZER_SETUP_FREE << std::endl;
        std::cout << "   - " << TWEEZER_SETUP_LINKFIX << std::endl;
        std::cout << "   - " << TWEEZER_SETUP_TORSIONAL_TRAP << std::endl;
        std::cout << "   - " << TWEEZER_SETUP_TORQUE << std::endl;
        std::exit(0);
    }





    std::cout << "Tweezer setup: " << tweezer_setup << std::endl << std::endl;


    /*
        Slow Winding
    */
    if (input->contains_multi("slow_wind")) {
        MultiLineInstance * ML = input->get_multiline("slow_wind")->get_instance(0);
        if (ML->contains_singleline("dLK_step") && ML->contains_singleline("steps") ) {
            slow_wind = true;
            slow_wind_dLK_step           = ML->get_single_val<double>("dLK_step",slow_wind_dLK_step);
            slow_wind_steps_per_dLK_step = ML->get_single_val<long long>("steps",slow_wind_steps_per_dLK_step);
            slow_wind_dump               = ML->get_single_val<bool>("dump",slow_wind_dump);

            std::cout << "Imposing dLK by slow winding: "  << std::endl;
            std::cout << " dLK step  =  " << slow_wind_dLK_step << std::endl;
            std::cout << " MC steps  =  " << slow_wind_steps_per_dLK_step << std::endl;
        }
    }

    check_link       = true;
    check_link       = InputChoice_get_single<bool> ("check_link",input,argv,check_link);
    std::cout << "check_link =  " << check_link << std::endl << std::endl;


    /*
        Chain Initialization
    */
    std::cout << "Initializing Chain .. " << std::endl;
    if (chain_initialized) {
        delete chain;
    }
    Chain * nchain = new Chain(IDB_fn);
    chain = nchain;
    chain_initialized = true;

    chain->set_T(temp);
    chain->set_force(force,fdir);
    chain->set_T0_subtract(T0_subtract);
//    chain->set_helical_repeat_length(hel_rep_len);

//    chain->fix_termini_position();
    chain->fix_termini_orientation();

    if (restart_file=="") {
        chain->gen_linear(num_bp, seq, sigma,fdir);
    }
    else {
        chain->gen_conf_from_restart(restart_file,restart_snapshot,"linear");
    }

    std::cout << std::endl;
    if (tweezer_setup==TWEEZER_SETUP_FREE) {
        std::cout << "No Link Constraint imposed .. " << std::endl;
//        chain->fix_link(false);
        chain->allow_free_endbead_rotation();
    }
    if (tweezer_setup==TWEEZER_SETUP_LINKFIX) {
        std::cout << "Fixed Link Constraint imposed .. " << std::endl;
        chain->fix_link(true);
    }
    if (tweezer_setup==TWEEZER_SETUP_TORSIONAL_TRAP) {
//        chain->fix_link(false);
        std::cout << "Inposing Torsional Trap .. " << std::endl;
        chain->impose_torsional_trap(chain->get_dLK());
    }
    if (tweezer_setup==TWEEZER_SETUP_TORQUE) {
        std::cout << "Inposing Torque Constraint to tau = " << torque << std::endl;
//        chain->fix_link(false);
        chain->impose_torque(torque);
    }




    /*
        MCStep Initialization
    */
    for (unsigned i=0;i<MCSteps.size();i++) {
        delete MCSteps[i];
    }
    MCSteps.clear();

    if (dynamics) {
        std::cout << std::endl << "Initializing MC Moves" << std::endl;
//        MCS_PivCon* pivcon    = new MCS_PivCon(chain,seedseq,2,num_bp/10);
        MCS_PivCon* pivcon    = new MCS_PivCon(chain,seedseq,2,15,larger(2,num_bp-15),num_bp);

        MCSteps.push_back(pivcon);
        std::cout << " Added PivCon to MC Moves .. " << std::endl;
        MCS_CSrot*  rot       = new MCS_CSrot(chain,seedseq,2,larger(2,10));
        MCSteps.push_back(rot);
        MCSteps.push_back(rot);
        MCSteps.push_back(rot);
        MCSteps.push_back(rot);
        MCSteps.push_back(rot);
        MCSteps.push_back(rot);
        MCSteps.push_back(rot);
        MCSteps.push_back(rot);
        MCSteps.push_back(rot);
        MCSteps.push_back(rot);
        std::cout << " Added CSrot to MC Moves .. " << std::endl;
    }
    else {
        std::cout << std::endl << "Initializing MC Moves" << std::endl;
        MCS_PivCon* pivcon    = new MCS_PivCon(chain,seedseq,2,num_bp/10);
        MCSteps.push_back(pivcon);
        std::cout << " Added PivCon to MC Moves .. " << std::endl;
        MCS_CSrot*  rot       = new MCS_CSrot(chain,seedseq,2,larger(2,num_bp/4));
        MCSteps.push_back(rot);
        std::cout << " Added CSrot to MC Moves .. " << std::endl;
    }

    if (use_cluster_twist) {
        MCS_ClusterTwist* cltwist = new MCS_ClusterTwist(chain,seedseq,2,larger(2,num_bp/4));
        MCSteps.push_back(cltwist);
        std::cout << " Added ClusterTwist to MC Moves .. " << std::endl;
    }

    if (tweezer_setup!=TWEEZER_SETUP_LINKFIX) {
        bool use_torpiv = InputChoice_get_single<bool> ("use_torpiv",input,argv,false);
        geninfile.add_entry(GENINFILE_SIMSETUP,"use_torpiv",use_torpiv);
        if (use_torpiv) {
            MCS_TorquePiv* torpiv = new MCS_TorquePiv(chain,seedseq);
            MCSteps.push_back(torpiv);
            std::cout << " Added torquePiv to MC Moves .. " << std::endl;
        }
        else {
            MCS_TailTwist* tailtw = new MCS_TailTwist(chain,seedseq);
            MCSteps.push_back(tailtw);
            std::cout << " Added TailTwist to MC Moves .. " << std::endl;
        }
    }

    if (num_twist>0) {
        MCS_Twist* twi = new MCS_Twist(chain,seedseq);
        for (unsigned i=0;i<num_twist;i++) {
            MCSteps.push_back(twi);
        }
        std::cout << " Added " << num_twist << " Twist-Moves to MC Moves .. " << std::endl;
    }

    ////////////////////////////////////////////
    // Kink move
    bool use_kinkxy = InputChoice_get_single<bool> ("use_kinkxy",input,argv,false);
    if (use_kinkxy) {
        // check if a BPS has energy eval type kinkxy
        bool kinkxy_included=false;
        for (unsigned i=0;i<chain->get_BPS()->size();i++) {
            (*chain->get_BPS())[i]->get_method_diag();
            if (*(*chain->get_BPS())[i]->get_method_diag() == "kinkxy") {
                kinkxy_included=true;
                break;
            }
        }
        if (kinkxy_included) {
            MCS_Kinkxy* mcs_kinkxy = new MCS_Kinkxy(chain,seedseq);
            MCSteps.push_back(mcs_kinkxy);
            std::cout << " Added Kinkxy to MC Moves .. " << std::endl;
        }
    }

    return true;
}


bool PolyMC::tweezer_slow_wind(double dLK_from, double dLK_to, double dLK_step, long long steps_per_dLK_step) {

    std::vector<MCStep*> wind_MCSteps;
    MCS_PivCon* pivcon    = new MCS_PivCon(chain,seedseq,2,num_bp/10);
    wind_MCSteps.push_back(pivcon);
    MCS_CSrot*  rot       = new MCS_CSrot(chain,seedseq,2,larger(2,num_bp/4));
    wind_MCSteps.push_back(rot);
    if (use_cluster_twist) {
        MCS_ClusterTwist* cltwist = new MCS_ClusterTwist(chain,seedseq,2,larger(2,num_bp/4));
        wind_MCSteps.push_back(cltwist);
    }

    if (EV_active) {
        for (unsigned i=0;i<wind_MCSteps.size();i++) {
            wind_MCSteps[i]->set_excluded_volume(EV);
        }
    }

    int sgn_of_change = sgn(dLK_to-dLK_from);
    if (sgn(dLK_step) != sgn_of_change) {
        dLK_step = -dLK_step;
    }

    double set_dLK = dLK_from;
    chain->set_Delta_Lk(set_dLK);
    EV->set_current_as_backup();

    if (check_link) {
        set_link_backup();
    }

    while (set_dLK*sgn_of_change < dLK_to*sgn_of_change ) {
        std::cout << "dLK set to " << set_dLK << std::endl;
        run(steps_per_dLK_step,wind_MCSteps,"winding chain - target: "+to_string_with_precision(dLK_to,2)+" turns\n",slow_wind_dump);
        set_dLK += dLK_step;

        if (set_dLK*sgn_of_change > dLK_to*sgn_of_change) {
            set_dLK = dLK_to;
        }
        chain->set_Delta_Lk(set_dLK);
        EV->set_current_as_backup();
    }
    std::cout << "dLK set to " << set_dLK << std::endl;
    run(steps_per_dLK_step,wind_MCSteps,"winding chain",slow_wind_dump);

    double Wr,Tw;
    Wr  = chain->cal_langowski_writhe_1a();
    Tw  = chain->cal_twist(0,num_bp);
    std::cout << "  LK = " << Tw+Wr << " ( " << chain->get_dLK() << ")" << std::endl;

    for (unsigned i=0;i<wind_MCSteps.size();i++) {
        delete wind_MCSteps[i];
    }
    return true;
}
