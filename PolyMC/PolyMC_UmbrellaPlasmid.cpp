#include "../PolyMC.h"

bool PolyMC::init_umbrellaplasmid() {

    std::cout << std::endl;
    std::cout << "######################################" << std::endl;
    std::cout << "## Loading Umbrella Plasmid Protocol #" << std::endl;
    std::cout << "######################################" << std::endl;
    std::cout << std::endl;

    closed                  = false;
    terminus_fixed_tangents = true;
    terminus_fixed_triads   = false;
    EV_repulsion_plane      = false;

    /*
        Umbrellaplasmid Settings
    */
    umbrellaplasmid_setup_id = -1;
    umbrellaplasmid_setup    = InputChoice_get_single<std::string> ("umbrellaplasmid_setup",input,argv,umbrellaplasmid_setup);
    geninfile.add_entry(GENINFILE_SIMSETUP,"umbrellaplasmid_setup",umbrellaplasmid_setup);
    if (umbrellaplasmid_setup==PLEC_SETUP_FREE) {
        umbrellaplasmid_setup_id = PLEC_SETUP_ID_FREE;
    }
    if (umbrellaplasmid_setup==PLEC_SETUP_LINKFIX) {
        umbrellaplasmid_setup_id = PLEC_SETUP_ID_LINKFIX;
        terminus_fixed_triads   = true;
    }
    if (umbrellaplasmid_setup==PLEC_SETUP_TORSIONAL_TRAP) {
        umbrellaplasmid_setup_id = PLEC_SETUP_ID_TORSIONAL_TRAP;
    }
    if (umbrellaplasmid_setup==PLEC_SETUP_TORQUE) {
        umbrellaplasmid_setup_id = PLEC_SETUP_ID_TORQUE;
        sigma = 0;
        std::cout << "Warning: Umbrellaplasmid setup Torque requires the initial supercoiling" << std::endl;
        std::cout << "density to be set to zero. (sigma was set to 0)" << std::endl;
    }
    if (umbrellaplasmid_setup_id<0) {
        std::cout << "Error: '" << umbrellaplasmid_setup << "' is not a valid umbrellaplasmid setup!" << std::endl;
        std::cout << " Valid setups:" << std::endl;
        std::cout << "   - " << PLEC_SETUP_FREE << std::endl;
        std::cout << "   - " << PLEC_SETUP_LINKFIX << std::endl;
        std::cout << "   - " << PLEC_SETUP_TORSIONAL_TRAP << std::endl;
        std::cout << "   - " << PLEC_SETUP_TORQUE << std::endl;
        std::exit(0);
    }
    std::cout << "Umbrellaplasmid setup: " << umbrellaplasmid_setup << std::endl << std::endl;

    /*
        Slow Winding
        TODO: This could be called for all modes. Within the mode init slow wind could get deactivated if
              not supported by the given setup.
    */
    if (input->contains_multi("slow_wind")) {
        std::cout << "Error: mode umbrellaplasid does not support slow winding" << std::endl;
        std::exit(0);
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
    chain->set_T0_subtract(T0_subtract);
    chain->set_helical_repeat_length(hel_rep_len);

    chain->fix_termini_position();
    chain->fix_termini_orientation();


    if (restart_file=="") {
        chain->gen_circular(num_bp, sigma, seq,false);
        chain->set_pseudo_closed_topology();
    }
    else {
        std::cout << "Error: Restart function not yet implemented for mode umbrellaplasmid." << std::endl;
        std::exit(0);
    }

    //                                       LK_0                              * sigma
    umbrellaplasmid_dLK   = (num_bp-1)*chain->get_disc_len()/hel_rep_len * sigma;
    chain->set_Delta_Lk(umbrellaplasmid_dLK);


    std::cout << std::endl;
    if (umbrellaplasmid_setup==PLEC_SETUP_FREE) {
        std::cout << "No Link Constraint imposed .. " << std::endl;
//        chain->fix_link(false);
        chain->allow_free_endbead_rotation();
    }
    if (umbrellaplasmid_setup==PLEC_SETUP_LINKFIX) {
        std::cout << "Fixed Link Constraint imposed .. " << std::endl;
        chain->fix_link(true);
    }
    if (umbrellaplasmid_setup==PLEC_SETUP_TORSIONAL_TRAP) {
//        chain->fix_link(false);
        std::cout << "Inposing Torsional Trap .. " << std::endl;
        umbrellaplasmid_trap_stiffness = InputChoice_get_single<double> ("umbrellaplasmid_trap_stiffness",input,argv,umbrellaplasmid_trap_stiffness);
        if (umbrellaplasmid_trap_stiffness <= 0) {
            std::cout << "Error: Torsional trap for Umbrellaplasmid mode requires the torsional trapstiffness to be set." << std::endl;
            std::cout << "  Specify trapstiffness with flag 'umbrellaplasmid_trap_stiffness'" << std::endl;
            std::exit(0);
        }
        chain->impose_torsional_trap(umbrellaplasmid_dLK,umbrellaplasmid_trap_stiffness);
    }
    if (umbrellaplasmid_setup==PLEC_SETUP_TORQUE) {
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
    MCS_CSrot*  rot       = new MCS_CSrot(chain,seedseq,2,larger(2,num_bp/4),0,num_bp-2);
    MCSteps.push_back(rot);
    MCS_CSrot*  rot2      = new MCS_CSrot(chain,seedseq,2,num_bp-1,0,num_bp-2);
    MCSteps.push_back(rot2);
    std::cout << " Added CSrot to MC Moves .. " << std::endl;
    if (use_cluster_twist) {
        MCS_ClusterTwist* cltwist = new MCS_ClusterTwist(chain,seedseq,2,larger(2,num_bp/4),1,num_bp-2);
        MCSteps.push_back(cltwist);
        std::cout << " Added ClusterTwist to MC Moves .. " << std::endl;
    }

    if (plec_setup!=PLEC_SETUP_LINKFIX) {
        MCS_TailTwist* tailtw = new MCS_TailTwist(chain,seedseq,1);
        MCSteps.push_back(tailtw);
        std::cout << " Added TailTwist to MC Moves .. " << std::endl;
    }

    if (num_twist>0) {
        MCS_Twist* twi = new MCS_Twist(chain,seedseq,1,num_bp-2);
        for (unsigned i=0;i<num_twist;i++) {
            MCSteps.push_back(twi);
        }
        std::cout << " Added " << num_twist << " Twist-Moves to MC Moves .. " << std::endl;
    }
    return true;
}
