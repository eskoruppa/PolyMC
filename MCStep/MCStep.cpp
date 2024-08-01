#include "MCStep.h"

MCStep::MCStep(Chain * ch,const std::vector<long long int> & seedseq) :
// Set State Variables
chain(ch),
seedseq(seedseq),
BPS(*ch->get_BPS()),
pos(ch->get_bp_pos()),
triads(ch->get_triads()),
states(ch->get_states()),
num_bp(ch->get_num_bp()),
num_bps(ch->get_num_bps()),
disc_len(ch->get_disc_len()),
count_step(0),
count_accept(0),
ev_active(false),
requires_EV_check(true),
es_active(false)
{
    // set seeds for random number generators
    std::seed_seq seed(seedseq.begin(), seedseq.end());
//    std::seed_seq seed{rd(), rd(), rd(), rd(), rd(), rd(), rd(), rd()};
    gen.seed(seed);
}


MCStep::~MCStep() {

}

bool MCStep::MC() {
    count_step++;

    if (gen_full_trial_conf) {

        double beta_Delta_E;
        bool   accepted = true;

        set_trial_backup();
        beta_Delta_E = gen_trial_conf();

        if (beta_Delta_E > 1e10) {
            accepted = false;
        }

        if (accepted && constraints_active && requires_constraint_check) {
            /*
                Check Constraints
            */
            accepted = check_constraints();
        }
        if (accepted && ev_active) {
            /*
                Check Excluded Volume
            */
            if (requires_EV_check) {
                accepted = EV->check(&moved_intervals);
                /*
                    currently EV automatically sets the backup if EV check is accepted
                */
            }
            else {
                /*
                TODO:   Build in option to selectively set backup to save time
                        by not copying the entire configuration
                */
                EV->set_current_as_backup();
            }
        }

        if (accepted) {
            if (pair_interactions_active) {
                /*
                    Calculate Pair Interactions
                */
                double pair_beta_delta_E;
                pair_beta_delta_E = chain->get_beta() * unbound->Eval_Delta_Energy(pos,&trial_backup_pos,moved_intervals);
                beta_Delta_E     += pair_beta_delta_E;
            }

            if (es_active) {
                /*
                    Calculate Electrostatics
                */
                beta_Delta_E += ES->cal_beta_dE(&moved_intervals);
            }

            /*
                Metropolis Step
            */
            if (std::exp(-beta_Delta_E) <= uniformdist(gen)) {
                accepted = false;

                revert_to_trial_backup();
                if (ev_active) {
                    /*
                        currently EV automatically sets the backup if EV check is accepted
                        ->  hence if move gets rejected on Metropolis step, we have to revert to
                            previous backup
                        TODO:
                            THIS SHOULD BE HANDLED OUTSIDE THE FUNCTION TO AVOID UNNECESSARY COPYING!
                    */
                    EV->set_current_as_backup();
                }
                if (es_active) {
                    ES->revert_to_backup();
                }
            }
            else {
                if (es_active) {
                    ES->set_current_as_backup();
                }
                count_accept++;
            }
        }
        else {
            revert_to_trial_backup();
        }

        for (unsigned i=0;i<changed_bps.n_elem;i++) {
            BPS[changed_bps(i)]->set_move(accepted);
        }

        if (!chain->check_energy_consistency())  {
            std::cout << "MCStep energy inconsistent" << std::endl;
            std::cout << "accepted " << accepted << std::endl;
            std::cout << "stepname " << move_name << std::endl;


            std::cout << "diff pos = " << arma::accu(*pos - trial_backup_pos) << std::endl;
            std::cout << "diff trd = " << arma::accu(*triads - trial_backup_triads) << std::endl;

            std::exit(0);
        }
        return accepted;
    }
    else {
        bool accepted = MC_move();
        if (accepted) {
            count_accept++;

            if (constraints_active && requires_constraint_check) {
                /*
                    Check Constraints
                */
                accepted = check_constraints();
            }
            if (!accepted) {
                count_accept--;
            }
            else if (ev_active) {
                /*
                    Check Excluded Volume
                */
                if (requires_EV_check) {
                    accepted = EV->check(&moved_intervals);
                    if (!accepted) {
                        count_accept--;
                    }
                }
                else {
                    /*
                    TODO:   Build in option to selectively set backup to save time
                            by not copying the entire configuration
                    */
                    EV->set_current_as_backup();
                }
            }
        }
        for (unsigned i=0;i<changed_bps.n_elem;i++) {
            BPS[changed_bps(i)]->set_move(accepted);
        }

        #ifdef MCSTEP_CONSISTENCYCHECK
        if (!check_consistency_violation()) {
            std::exit(0);
        }
        #endif

        return accepted;
    }
}

bool MCStep::MC_move() {
    /*
        Here the corpus of the monte carlo step should go
    */
    std::cout << "Error: The base class MCStep should never be called for MC moves." << std::endl;
    std::exit(0);
    return true;
}


void MCStep::set_excluded_volume(ExVol* EVol) {
    ev_active = true;
    EV        = EVol;
//    if (requires_EV_check) {
//
//    }
//    else {
//        cout << move_name << " does not require excluded volume interactions or is incompatible with excluded volumes." << endl;
//    }
}

void MCStep::set_constraints(const std::vector<Constraint*> & constr) {
    if (constr.size() > 0) {
        constraints = constr;
        constraints_active = true;
    }
}

void MCStep::remove_constraints() {
    constraints.empty();
    constraints_active = false;
}


std::vector<arma::ivec>* MCStep::get_moved_intervals() {
    return &moved_intervals;
}

/*--------------------------------------------------------------------*/
/*--------------------------------------------------------------------*/
/*
    Electro Statics methods
*/
/*--------------------------------------------------------------------*/


void MCStep::set_electrostatics(ElStat * elstat) {
    es_active           = true;
    ES                  = elstat;
    if (requires_EV_check){
        gen_full_trial_conf = true;
    }
}

/*--------------------------------------------------------------------*/
/*--------------------------------------------------------------------*/





bool MCStep::suitable() {
    // force_active
    if (chain->force_active() && !suitablility_key[0]) {
        std::cout << "force conflict" << std::endl;
        std::cout << suitablility_key[0] << std::endl;
        return false;
    }
    if (chain->const_torque_active() && !suitablility_key[1]) {
        return false;
    }
    if (chain->topology_closed() && !suitablility_key[2]) {
        return false;
    }
    if (chain->link_fixed() && !suitablility_key[3]) {
        return false;
    }
    if (chain->fixed_termini() && !suitablility_key[4]) {
        return false;
    }
    if (chain->fixed_termini_radial() && !suitablility_key[5]) {
        return false;
    }
    if (chain->fixed_first_orientation() && !suitablility_key[6]) {
        return false;
    }
    if (chain->fixed_last_orientation() && !suitablility_key[7]) {
        return false;
    }
    return additional_criteria();
}

bool MCStep::additional_criteria() {
    /*
        Additional criteria beyond the standard criteria to verify if this monte carlo
        Move is applicable under the given conditions of the system. Should be overwritten
        by MC moves that inherit MCStep.
    */
    return true;
}

void MCStep::reset_count() {
    count_ev_accept = 0;
    count_step      = 0;
    count_accept    = 0;
}


void MCStep::update_settings() {}


long long int MCStep::get_steps() {
    return count_step;
}
long long int MCStep::get_accepts() {
    return count_accept;
}
long long int MCStep::get_ev_accepts() {
    return count_ev_accept;
}
std::string   MCStep::get_move_name() {
    return move_name;
}
double MCStep::acceptance_rate() {
    if (count_accept>0) {
        return (double)count_accept/(double)count_step;
    }
    else {
        return 0;
    }
}
void MCStep::print_acceptance_rate() {
    std::string acrte = " Acceptance Rate ";
    std::string outstr = "  " + move_name;
    int diff = MSC_ACCEPTANCE_PRINT_STRLEN-outstr.length()-acrte.length();
    for (int i=0;i<diff;i++) {
        outstr += " ";
    }
    outstr += acrte;
    std::cout << outstr << " = " << acceptance_rate() << std::endl;
}


/*--------------------------------------------------------------------*/
/*--------------------------------------------------------------------*/
/*
    Check Constraints
*/

bool MCStep::check_constraints() {
    for (unsigned cstr=0;cstr<constraints.size();cstr++){
        if (!constraints[cstr]->check(&moved_intervals)) {
            return false;
        }
    }
    return true;
}


/*--------------------------------------------------------------------*/



/*--------------------------------------------------------------------*/
/*--------------------------------------------------------------------*/
/*
    Pair Interactions
*/
/*--------------------------------------------------------------------*/



void MCStep::set_unbound(Unbound * unb) {
    unbound = unb;
    if (requires_pair_check) {
        pair_interactions_active = true;
        gen_full_trial_conf      = true;
        std::cout << "Set pair interactions for " << move_name << std::endl;
    }
    else {
        pair_interactions_active = false;
    }
//    std::cout << " paircheck required: " << move_name << " - " << requires_pair_check << std::endl;
}

void MCStep::set_unbound_active(bool active) {
    pair_interactions_active = active;
}


/*--------------------------------------------------------------------*/
/*--------------------------------------------------------------------*/
/*
    Variables and Methods for Full Configuration generatation (used for Pair interactions)
*/
/*--------------------------------------------------------------------*/

double MCStep::gen_trial_conf() {
    /*
        Here the corpus of the monte carlo step should go
    */
    std::cout << "Error: The base class MCStep should never be called to generate trial configrations." << std::endl;
    std::exit(0);
    return 0;
}

void MCStep::set_trial_backup() {
    trial_backup_pos    = *pos;
    trial_backup_triads = *triads;
}
void MCStep::revert_to_trial_backup() {
    *pos    = trial_backup_pos;
    *triads = trial_backup_triads;
}

/*--------------------------------------------------------------------*/
/*--------------------------------------------------------------------*/
/*
    Consistency checks for debugging
*/
/*--------------------------------------------------------------------*/

bool MCStep::check_consistency_violation(bool print_consistent) {
    if (!chain->check_energy_consistency()) {
        std::cout << move_name << " -> Energy inconsistent!" << std::endl;
        return false;
    }
    else if (print_consistent) {
        std::cout << move_name << " -> Energy consistent!" << std::endl;
    }
    if (!chain->config_consistent()) {
        std::cout << move_name << "-> chain positions inconsistent!" << std::endl;
        return false;
    }
    else if (print_consistent) {
        std::cout << move_name << "-> chain positions consistent!" << std::endl;
    }
    return true;
}




