#include "MCS_Unwinder.h"


MCS_Unwinder::MCS_Unwinder(Chain * ch,const std::vector<long long int> & seedseq, double sigma_unwound, int steps_unwinding, int steps_winding, int step_equi, int n_split_winding)
: MCStep(ch,seedseq),
sigma_unwound(sigma_unwound),
steps_unwinding(steps_unwinding),
steps_winding(steps_winding),
steps_winding_per_split(steps_winding/n_split_winding),
step_equi(step_equi),
n_split_winding(n_split_winding)
{
    move_name = "MCS_Unwinder";

    dLK_unwound = chain->sigma2dLk(sigma_unwound);
    dLK_wound   = chain->get_dLK();

    double dLK_increment   = (double)((dLK_wound-dLK_unwound))/n_split_winding;
    double dLK_vals        = dLK_unwound;
    for (int i=0;i<n_split_winding;i++) {
        dLK_vals += dLK_increment;
        winding_dLK.push_back(dLK_vals);
    }

    MCS_CSrot * rot1 = new MCS_CSrot(chain,seedseq,2,larger(2,num_bp/2));
    MCS_CSrot * rot2 = new MCS_CSrot(chain,seedseq,2,larger(2,num_bp/2));
    MCS_CSrot * rot3 = new MCS_CSrot(chain,seedseq,2,larger(2,num_bp/2));
    MCS_CSrot * rot4 = new MCS_CSrot(chain,seedseq,2,larger(2,num_bp/2));
    MCS_CSrot * rot5 = new MCS_CSrot(chain,seedseq,2,larger(2,num_bp/2));

    MCSteps.push_back(rot1);
    MCSteps.push_back(rot2);
    MCSteps.push_back(rot3);
    MCSteps.push_back(rot4);
    MCSteps.push_back(rot5);

    link_check_density = disc_len/MCSUW_CHECK_LINK_SEG_SIZE;
    if (link_check_density > 1) link_check_density = 1;

    backup_bp_pos = *chain->get_bp_pos();
    backup_triads = *chain->get_triads();

    suitablility_key[MCS_FORCE_ACTIVE]              = false;  // force_active
    suitablility_key[MCS_TORQUE_ACTIVE]             = false;  // torque_active
    suitablility_key[MCS_TOPOLOGY_CLOSED]           = true;  // topology_closed
    suitablility_key[MCS_FIXED_LINK]                = true;  // link_fixed
    suitablility_key[MCS_FIXED_TERMINI]             = false;  // fixed_termini
    suitablility_key[MCS_FIXED_TERMINI_RADIAL]      = false;  // fixed_termini_radial
    suitablility_key[MCS_FIXED_FIRST_ORIENTATION]   = false;  // fixed_first_orientation
    suitablility_key[MCS_FIXED_LAST_ORIENTATION]    = false;  // fixed_last_orientation
}

MCS_Unwinder::~MCS_Unwinder() {
}

void MCS_Unwinder::set_excluded_volume(ExVol* EV) {
    this->EV = EV;
    for (unsigned i=0;i<MCSteps.size();i++) {
        MCSteps[i]->set_excluded_volume(EV);
    }
}

void MCS_Unwinder::update_settings() {
    for (unsigned i=0;i<MCSteps.size();i++) {
        MCSteps[i]->update_settings();
    }
}

bool MCS_Unwinder::MC() {
    bool accepted = MC_move();
    if (accepted) {
        count_accept++;
    }
    return accepted;
}


bool MCS_Unwinder::MC_move() {

    double E_init,E_final,E_diff;
    bool   accepted;

    /*
        Store initial energy and make a copy of the configuration
    */
    E_init = chain->extract_energy();
    backup_bp_pos = *chain->get_bp_pos();
    backup_triads = *chain->get_triads();

    /*
        Unwinding
    */
    std::cout << "Unwinding..." << std::endl;
    std::cout << "dLK set to " << dLK_unwound << std::endl;
    chain->set_Delta_Lk(dLK_unwound);
    EV->set_backup_conf(chain->get_bp_pos(),chain->get_triads());
    run(steps_unwinding);

    /*
        Winding
    */
    std::cout << "Winding..." << std::endl;
    for (int w=0;w<n_split_winding;w++) {
        std::cout << "dLK set to " << winding_dLK[w] << std::endl;
        chain->set_Delta_Lk(winding_dLK[w]);
        EV->set_backup_conf(chain->get_bp_pos(),chain->get_triads());
        run(steps_winding_per_split);
    }

    /*
        Equilibrate
    */
    std::cout << "Equilibrate..." << std::endl;
    chain->set_Delta_Lk(dLK_wound);
    EV->set_backup_conf(chain->get_bp_pos(),chain->get_triads());
    run(step_equi);

    /*
        Metropolis
    */

    E_final = chain->extract_energy();

//    E_diff = E_final-E_init;
//    E_diff = (E_final-E_init)/std::sqrt(num_bps);
    E_diff = (E_final-E_init)/num_bps;

    double weight = exp(-E_diff);

    std::cout << "Energy Difference         = " << E_final-E_init << std::endl;
    std::cout << "Considered Difference     = " << E_diff << std::endl;
    std::cout << "Boltzmann Weight          = " << weight << std::endl;



    accepted = weight > uniformdist(gen);

    if (!accepted) {
        chain->set_config(&backup_bp_pos,&backup_triads,true);
        EV->set_backup_conf(&backup_bp_pos,&backup_triads);
    }
    return accepted;
}



void MCS_Unwinder::run(long long int steps) {
    for (long long int step=0;step<=steps;step++) {
        if (step%MCSUW_CHECK_LINK_EVERY==0) {
            if (!chain->check_link_conservation()) {
                if (MCSUW_KILL_ON_DLK_VIOLATION) {
                    std::cout << "Conservation of Linking Number violated!" << std::endl;
                    std::exit(0);
                }
            }
        }
        for (unsigned i=0;i<MCSteps.size();i++) {
            MCSteps[i]->MC();
        }
    }
}


