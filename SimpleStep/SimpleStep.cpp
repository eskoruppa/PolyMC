#include "MCStep.h"

MCStep::MCStep(Chain * ch) :
// Set State Variables
chain(ch),
BPS(*ch->get_BPS()),
pos(ch->get_bp_pos()),
triads(ch->get_triads()),
num_bp(ch->get_num_bp()),
num_bps(ch->get_num_bps()),
disc_len(ch->get_disc_len()),
count_step(0),
count_accept(0),
ev_active(false),
requires_EV_check(true)
{
    // set seeds for random number generators
    std::seed_seq seed{rd(), rd(), rd(), rd(), rd(), rd(), rd(), rd()};
    gen.seed(seed);
}


MCStep::~MCStep() {

}

bool MCStep::MC() {
    count_step++;
    bool accepted = MC_move();
    if (accepted) {
        count_accept++;
        if (ev_active) {
            /*
                Check Excluded Volume
            */
            accepted = EV->check(&moved_intervals);
            if (!accepted) {
                count_accept--;
                chain->recal_energy();
            }
        }
    }
    for (int i=0;i<changed_bps.n_elem;i++) {
        BPS[changed_bps(i)]->set_move(accepted);
    }
    return accepted;
}

bool MCStep::MC_move() {
    /*
        Here the corpus of the monte carlo step should go
    */
    return true;
}


void MCStep::set_excluded_volume(ExVol* EVol) {
    if (requires_EV_check) {
        ev_active = true;
        EV        = EVol;
    }
//    else {
//        cout << move_name << " does not require excluded volume interactions or is incompatible with excluded volumes." << endl;
//    }
}


vector<arma::ivec>* MCStep::get_moved_intervals() {
    return &moved_intervals;
}

bool MCStep::suitable() {
    // force_active
    if (chain->force_active() && !suitablility_key[0]) {
        cout << "force conflict" << endl;
        cout << suitablility_key[0] << endl;
        return false;
    }
    if (chain->torque_active() && !suitablility_key[1]) {
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
double MCStep::acceptance_rate() {
    if (count_accept>0) {
        return (double)count_accept/(double)count_step;
    }
    else {
        return 0;
    }
}
void MCStep::print_acceptance_rate() {
    cout << move_name << " Acceptance Rate = " << acceptance_rate() << endl;
}

