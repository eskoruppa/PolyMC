#include "../PolyMC.h"

bool PolyMC::init_open() {

    std::cout << std::endl;
    std::cout << "######################################" << std::endl;
    std::cout << "##### Loading Open Chain Protocol ####" << std::endl;
    std::cout << "######################################" << std::endl;
    std::cout << std::endl;

    closed                  = false;
    terminus_fixed_tangents = false;
    terminus_fixed_triads   = false;
    EV_repulsion_plane      = false;

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

    // set closure parameters
    chain->set_closure_distance_stiff(closure_distance_stiff,closure_distance_equi);
    chain->set_closure_angle_stiff(closure_angle_stiff,closure_angle_use_costheta);

    chain->set_T0_subtract(subtract_T0);
    chain->set_helical_repeat_length(hel_rep_len);
    chain->set_Lk0_from_static(Lk0_from_static);
//    chain->gen_linear(num_bp, seq, sigma,fdir);

    if (restart_file=="") {
        chain->gen_linear(num_bp, seq, sigma,fdir);
    }
    else {
//        std::cout << "LOADING RESTART" << std::endl;
//        std::exit(0);
        chain->gen_conf_from_restart(restart_file,restart_snapshot,"linear");
    }


    /*
        MCStep Initialization
    */
    for (unsigned i=0;i<MCSteps.size();i++) {
        delete MCSteps[i];
    }
    MCSteps.clear();

    std::cout << std::endl << "Initializing MC Moves" << std::endl;
    MCS_Pivot* piv = new MCS_Pivot(chain,seedseq);
    MCSteps.push_back(piv);
    std::cout << " Added Pivot to MC Moves .. " << std::endl;

    // // TEST CODE!
    // MCS_CSrot*  rot       = new MCS_CSrot(chain,seedseq,2,larger(2,num_bp/4));
    // MCSteps.push_back(rot);
    // std::cout << " Added CSrot to MC Moves .. " << std::endl;

    if (use_cluster_twist) {
        MCS_ClusterTwist* cltwist = new MCS_ClusterTwist(chain,seedseq,2,larger(2,num_bp/4));
        MCSteps.push_back(cltwist);
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
