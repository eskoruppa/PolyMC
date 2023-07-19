#include "../PolyMC.h"

bool PolyMC::init_lin2d() {

    std::cout << std::endl;
    std::cout << "######################################" << std::endl;
    std::cout << "##### Loading Linear2D Protocol ######" << std::endl;
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
    chain->set_T0_subtract(T0_subtract);
    chain->set_helical_repeat_length(hel_rep_len);
    chain->set_Lk0_from_static(Lk0_from_static);

//    chain->gen_linear(num_bp, seq, sigma,fdir);

    if (restart_file=="") {
        chain->gen_linear(num_bp, seq, sigma,fdir);
    }
    else {
        chain->gen_conf_from_restart(restart_file,restart_snapshot,"linear");
    }

    /*
        MCStep Initialization
    */
    for (unsigned i=0;i<MCSteps.size();i++) {
        delete MCSteps[i];
    }
    MCSteps.clear();

    use_slither2d = InputChoice_get_single<bool> ("use_slither2d",input,argv,use_slither2d);

    std::cout << std::endl << "Initializing MC Moves" << std::endl;
    if (use_slither2d) {
        MCS_Slither2d* sli2d = new MCS_Slither2d(chain,seedseq,2,larger(2,num_bp/4),{1,0,0},true);
        MCSteps.push_back(sli2d);
        std::cout << " Added MCS_Slither2d to MC Moves .. " << std::endl;
    }
    else {
        MCS_Pivot2d*   piv2d = new MCS_Pivot2d(chain,seedseq,{1,0,0});
        MCSteps.push_back(piv2d);
        std::cout << " Added MCS_Pivot2d   to MC Moves .. " << std::endl;
    }

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

bool PolyMC::lin2d_slither_equi_with_pivot() {

    std::vector<MCStep*> equi_MCSteps;
    MCS_Pivot2d*   piv2d = new MCS_Pivot2d(chain,seedseq,{1,0,0});
    equi_MCSteps.push_back(piv2d);

    if (EV_active) {
        for (unsigned i=0;i<equi_MCSteps.size();i++) {
            equi_MCSteps[i]->set_excluded_volume(EV);
        }
    }
    run(num_bp*100,equi_MCSteps,"Equil. slither2d");
    return true;
}

