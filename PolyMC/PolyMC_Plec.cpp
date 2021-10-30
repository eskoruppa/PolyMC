#include "../PolyMC.h"

bool PolyMC::init_plec() {

    std::cout << std::endl;
    std::cout << "######################################" << std::endl;
    std::cout << "#### Loading Plectoneme Protocol #####" << std::endl;
    std::cout << "######################################" << std::endl;
    std::cout << std::endl;

    closed                  = false;
    terminus_fixed_tangents = true;
    terminus_fixed_triads   = false;
    EV_repulsion_plane      = false;


    /*
        Plec Settings
    */
    plec_setup_id = -1;
    plec_setup        = InputChoice_get_single<std::string> ("plec_setup",input,argv,plec_setup);
    geninfile.add_entry(GENINFILE_SIMSETUP,"plec_setup",plec_setup);
    if (plec_setup==PLEC_SETUP_FREE) {
        plec_setup_id = PLEC_SETUP_ID_FREE;
    }
    if (plec_setup==PLEC_SETUP_LINKFIX) {
        plec_setup_id = PLEC_SETUP_ID_LINKFIX;
        terminus_fixed_triads   = true;
    }
    if (plec_setup==PLEC_SETUP_TORSIONAL_TRAP) {
        plec_setup_id = PLEC_SETUP_ID_TORSIONAL_TRAP;
    }
    if (plec_setup==PLEC_SETUP_TORQUE) {
        plec_setup_id = PLEC_SETUP_ID_TORQUE;
        sigma = 0;
        std::cout << "Warning: Plec setup Torque requires the initial supercoiling" << std::endl;
        std::cout << "density to be set to zero. (sigma was set to 0)" << std::endl;
    }
    if (plec_setup_id<0) {
        std::cout << "Error: '" << plec_setup << "' is not a valid plec setup!" << std::endl;
        std::cout << " Valid setups:" << std::endl;
        std::cout << "   - " << PLEC_SETUP_FREE << std::endl;
        std::cout << "   - " << PLEC_SETUP_LINKFIX << std::endl;
        std::cout << "   - " << PLEC_SETUP_TORSIONAL_TRAP << std::endl;
        std::cout << "   - " << PLEC_SETUP_TORQUE << std::endl;
        std::exit(0);
    }
    std::cout << "Plec setup: " << plec_setup << std::endl << std::endl;


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
//    chain->set_force(force,fdir);
    chain->set_T0_subtract(T0_subtract);
    chain->set_helical_repeat_length(hel_rep_len);

    chain->fix_termini_position();
    chain->fix_termini_orientation();


    if (restart_file=="") {
        plec_termini_dist = 2*EV_rad;
        plec_loop_frac    = 2./3;
        plec_termini_dist = InputChoice_get_single<double> ("plec_termini_dist",input,argv,plec_termini_dist);
        plec_loop_frac    = InputChoice_get_single<double> ("plec_loop_frac"   ,input,argv,plec_loop_frac   );

        arma::mat   pos;
        arma::cube  triads;

        std::cout << "Generating Plectoneme Configuration " << std::endl;
        std::cout << " Termini distance: " << plec_termini_dist << std::endl;
        std::cout << " Loop fraction:    " << plec_loop_frac << std::endl;

        plec_id_first      = gen_PFE_conf(&pos, &triads, num_bp, chain->get_disc_len(), plec_termini_dist, EV_rad, plec_loop_frac);
        num_bp             = pos.n_cols;
        plec_num_bp_actual = num_bp-plec_id_first;

        std::string addseq = "";
        for (int i=0;i<plec_id_first;i++) {
            addseq += seq[0];
        }
        seq = addseq+seq;
        chain->init_custom_open_conf(pos, triads, seq);
    }
    else {
        std::cout << "Error: Restart function not yet implemented for mode plec." << std::endl;
        std::exit(0);
//        chain->gen_conf_from_restart(restart_file,restart_snapshot,"linear");
    }

    //                                       LK_0                              * sigma
    plec_dLK            = plec_num_bp_actual*chain->get_disc_len()/hel_rep_len * sigma;
    chain->set_Delta_Lk(plec_dLK,plec_id_first,num_bp-1);

    std::cout << std::endl;
    if (plec_setup==PLEC_SETUP_FREE) {
        std::cout << "No Link Constraint imposed .. " << std::endl;
//        chain->fix_link(false);
        chain->allow_free_endbead_rotation();
    }
    if (plec_setup==PLEC_SETUP_LINKFIX) {
        std::cout << "Fixed Link Constraint imposed .. " << std::endl;
        chain->fix_link(true);
    }
    if (plec_setup==PLEC_SETUP_TORSIONAL_TRAP) {
//        chain->fix_link(false);
        std::cout << "Inposing Torsional Trap .. " << std::endl;
        plec_trap_stiffness = InputChoice_get_single<double> ("plec_trap_stiffness",input,argv,plec_trap_stiffness);
        if (plec_trap_stiffness <= 0) {
            std::cout << "Error: Torsional trap for Plectoneme mode requires the torsional trapstiffness to be set." << std::endl;
            std::cout << "  Specify trapstiffness with flag 'plec_trap_stiffness'" << std::endl;
            std::exit(0);
        }
        chain->impose_torsional_trap(plec_dLK,plec_trap_stiffness);
    }
    if (plec_setup==PLEC_SETUP_TORQUE) {
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

    std::cout << std::endl << "Initializing MC Moves" << std::endl;
    MCS_CSrot*  rot       = new MCS_CSrot(chain,seedseq,2,larger(2,num_bp/4),plec_id_first+1,num_bp-2);
    MCSteps.push_back(rot);
    MCS_CSrot*  rot2      = new MCS_CSrot(chain,seedseq,2,num_bp-3-plec_id_first,plec_id_first+1,num_bp-2);
    MCSteps.push_back(rot2);
    std::cout << " Added CSrot to MC Moves .. " << std::endl;
    if (use_cluster_twist) {
        MCS_ClusterTwist* cltwist = new MCS_ClusterTwist(chain,seedseq,2,larger(2,num_bp/4),plec_id_first+1,num_bp-2);
        MCSteps.push_back(cltwist);
        std::cout << " Added ClusterTwist to MC Moves .. " << std::endl;
    }

    if (plec_setup!=PLEC_SETUP_LINKFIX) {
        MCS_TailTwist* tailtw = new MCS_TailTwist(chain,seedseq,plec_id_first+1);
        MCSteps.push_back(tailtw);
        std::cout << " Added TailTwist to MC Moves .. " << std::endl;
    }

    if (num_twist>0) {
        MCS_Twist* twi = new MCS_Twist(chain,seedseq,plec_id_first+1,num_bp-2);
        for (unsigned i=0;i<num_twist;i++) {
            MCSteps.push_back(twi);
        }
        std::cout << " Added " << num_twist << " Twist-Moves to MC Moves .. " << std::endl;
    }
    return true;
}

bool PolyMC::plec_slow_wind(double dLK_from, double dLK_to, double dLK_step, long long steps_per_dLK_step) {

    std::cout << "################################################" << std::endl;
    std::cout << "################################################" << std::endl;
    std::cout << "############### Chain Winding ##################" << std::endl;
    std::cout << "Winding Chain from " << dLK_from << " turns to " << dLK_to << " turns in steps of " << dLK_step << " turns," << std::endl;
    std::cout << "running " << steps_per_dLK_step << " MC cycles per step." << std::endl;

    std::vector<MCStep*> wind_MCSteps;
    MCS_CSrot*  rot       = new MCS_CSrot(chain,seedseq,2,larger(2,num_bp/4),plec_id_first+1,num_bp-2);
    wind_MCSteps.push_back(rot);
    if (use_cluster_twist) {
        MCS_ClusterTwist* cltwist = new MCS_ClusterTwist(chain,seedseq,2,larger(2,num_bp/4),plec_id_first+1,num_bp-2);
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
    chain->set_Delta_Lk(set_dLK,plec_id_first,num_bp-1);
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
        chain->set_Delta_Lk(set_dLK,plec_id_first,num_bp-1);
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


