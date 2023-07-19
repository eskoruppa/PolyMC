#include "../PolyMC.h"

bool PolyMC::init_closed2d() {

    std::cout << std::endl;
    std::cout << "######################################" << std::endl;
    std::cout << "##### Loading Closed 2D Protocol #####" << std::endl;
    std::cout << "######################################" << std::endl;
    std::cout << std::endl;

    closed                  = true;
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

    chain->fix_link(true);
    chain->set_T(temp);
    chain->set_T0_subtract(T0_subtract);
    chain->set_helical_repeat_length(hel_rep_len);
    chain->set_Lk0_from_static(Lk0_from_static);
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
    MCS_Slither2d* sli2d = new MCS_Slither2d(chain,seedseq,2,larger(2,num_bp/4),{0,0,1});
    MCSteps.push_back(sli2d);
    std::cout << " Added Slither2d to MC Moves .. " << std::endl;

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
