#include "../Chain.h"

//////////////////////////////////////////////////////////
//////////// CALCULATE AND ACCESS ENERGIES ///////////////
//////////////////////////////////////////////////////////

/*
    These functions allow to calculate the new energy of
    a new configuration that is already stored in the
    triads
*/

void Chain::cal_energy_propose (int from, int to) {
    if (from>to && closed_topology) {
        cal_energy_propose(0,to);
        cal_energy_propose(from,num_bp-1);
        return;
    }
    for (int i=from;i<to;i++) {
        BPS[i]->propose_move(triads.slice(i),triads.slice(pmod(i+1,num_bp)));
    }
    BPS[to]->propose_move(triads.slice(to),triads.slice(pmod(to+1,num_bp)));
}

double Chain::cal_energy_eval (int from, int to) {
    double energy=0;
    if (from>to && closed_topology) {
        energy += cal_energy_eval(0,to);
        energy += cal_energy_eval(from,num_bp-1);
        return energy;
    }
    for (int i=from;i<=to;i++) {
        energy += BPS[i]->eval_new_energy();
    }
    return energy;
}

double Chain::cal_energy (int from, int to) {
    cal_energy_propose(from,to);
    return cal_energy_eval(from,to);
}

void   Chain::cal_energy_set_move(int from, int to,bool accept) {
    if (from>to && closed_topology) {
        cal_energy_set_move(0,to,accept);
        cal_energy_set_move(from,num_bp-1,accept);
        return;
    }
    for (int i=from;i<=to;i++) {
        BPS[i]->set_move(accept);
    }
}

void   Chain::set_new_energy() {
    cal_energy_propose(0,num_bps);
    cal_energy_eval(0,num_bps);
    cal_energy_set_move(0,num_bps,true);
}

//////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////
/*
    With those functions current energies in a given range
    can be read out.
*/

void  Chain::extract_energy_select(int from, int to) {
    if (from>to && closed_topology) {
        extract_energy_select(0,to);
        extract_energy_select(from,num_bp-1);
        return;
    }
    for (int i=from;i<=to;i++) {
        BPS[i]->energy_extract_select();
    }
}

double  Chain::extract_energy() {
    return extract_energy(0,num_bps-1);
}

double  Chain::extract_energy(int from, int to) {
    extract_energy_select(from,to);
    return __extract_energy(from,to);
}

double  Chain::__extract_energy(int from, int to) {
    double energy=0;
    if (from>to && closed_topology) {
        energy  = __extract_energy(0,to);
        energy += __extract_energy(from,num_bps-1);
        return energy;
    }
    for (int i=from;i<=to;i++) {
        energy += BPS[i]->energy_extract();
    }
    return energy;
}

//////////////////////////////////////////////////////////
////////////// Calculate True Energy /////////////////////
//////////////////////////////////////////////////////////

double  Chain::extract_true_energy() {
    return extract_energy(0,num_bps-1)/beta;
}

//////////////////////////////////////////////////////////
////////////// Calculate Full Energy /////////////////////
//////////////////////////////////////////////////////////

double Chain::extract_full_energy() {
    double E = extract_true_energy();
    E += extract_torque_betaenergy()/beta;
    E += extract_force_betaenergy()/beta;
    return E;
}


//////////////////////////////////////////////////////////
////////////// CHECK ENERGY CONSISTENCY //////////////////
//////////////////////////////////////////////////////////

bool Chain::check_energy_consistency() {
    double E_old, E_new;
    E_old = extract_energy();

    if (closed_topology) {
        for (unsigned bps=0;bps<num_bps-1;bps++) {
            BPS[bps]->propose_move(triads.slice(bps),triads.slice(bps+1));
        }
        BPS[num_bps-1]->propose_move(triads.slice(num_bps-1),triads.slice(0));
    }
    else {
        for (unsigned bps=0;bps<num_bps;bps++) {
            BPS[bps]->propose_move(triads.slice(bps),triads.slice(bps+1));
        }
    }
    for (unsigned bps=0;bps<num_bps;bps++) {
        BPS[bps]->eval_delta_energy();
    }
    for (unsigned bps=0;bps<num_bps;bps++) {
        BPS[bps]->set_move(true);
    }
    E_new = extract_energy();
    if (!(std::abs(E_new-E_old)<MAX_ENERGY_DISCREPANY)) {
        std::cout << "Energy difference = " << std::abs(E_new-E_old) << std::endl;
    }
    return (std::abs(E_new-E_old)<MAX_ENERGY_DISCREPANY);
}


double Chain::recal_energy() {

    if (closed_topology) {
        for (unsigned bps=0;bps<num_bps-1;bps++) {
            BPS[bps]->propose_move(triads.slice(bps),triads.slice(bps+1));
        }
        BPS[num_bps-1]->propose_move(triads.slice(num_bps-1),triads.slice(0));
    }
    else {
        for (unsigned bps=0;bps<num_bps;bps++) {
            BPS[bps]->propose_move(triads.slice(bps),triads.slice(bps+1));
        }
    }
    for (unsigned bps=0;bps<num_bps;bps++) {
        BPS[bps]->eval_delta_energy();
    }
    for (unsigned bps=0;bps<num_bps;bps++) {
        BPS[bps]->set_move(true);
    }
    return extract_energy();
}
