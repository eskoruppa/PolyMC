#include "../PolyMC.h"

bool PolyMC::init_plasmid() {

    std::cout << std::endl;
    std::cout << "######################################" << std::endl;
    std::cout << "###### Loading Plasmid Protocol ######" << std::endl;
    std::cout << "######################################" << std::endl;
    std::cout << std::endl;

    closed                  = true;
    terminus_fixed_tangents = false;
    terminus_fixed_triads   = false;
    EV_repulsion_plane      = false;


    /*
        Slow Winding
        TODO: This could be called for all modes. Within the mode init slow wind could get deactivated if
              not supported by the given setup.
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

    chain->fix_link(true);
    chain->set_T(temp);
    chain->set_T0_subtract(T0_subtract);
    chain->set_helical_repeat_length(hel_rep_len);
//    chain->gen_circular(num_bp, sigma, seq);

    if (restart_file=="") {
        chain->gen_circular(num_bp, sigma, seq);
    }
    else {
        chain->gen_conf_from_restart(restart_file,restart_snapshot,"circular");
    }

    /*
        MCStep Initialization
    */
    for (unsigned i=0;i<MCSteps.size();i++) {
        delete MCSteps[i];
    }
    MCSteps.clear();

    std::cout << std::endl << "Initializing MC Moves" << std::endl;

    if (InputChoice_get_single<int> ("dynamics",input,argv,0) > 0) {
        int dyn_dist = InputChoice_get_single<int> ("dynamics",input,argv,0);
        std::cout << " Dynamic moves set to " << dyn_dist << std::endl;
        MCS_CSrot* rot = new MCS_CSrot(chain,seedseq,2,larger(2,smaller(dyn_dist,num_bp/2)));
        MCSteps.push_back(rot);
    }
    else {
        MCS_CSrot* rot = new MCS_CSrot(chain,seedseq,2,larger(2,num_bp/2));
        MCSteps.push_back(rot);
    }
//    MCS_CSrot* rot = new MCS_CSrot(chain,2,100);
    std::cout << " Added CSrot to MC Moves .. " << std::endl;

    if (use_cluster_twist) {
        MCS_ClusterTwist* cltw = new MCS_ClusterTwist(chain,seedseq,2,larger(2,num_bp/2));
        MCSteps.push_back(cltw);
        std::cout << " Added ClusterTwist to MC Moves .. " << std::endl;
    }

    if (num_twist>0) {
        MCS_Twist* twi = new MCS_Twist(chain,seedseq);
        for (unsigned i=0;i<num_twist;i++) {
            MCSteps.push_back(twi);
        }
        std::cout << " Added " << num_twist << " Twist-Moves to MC Moves .. " << std::endl;
    }
    return true;
}
