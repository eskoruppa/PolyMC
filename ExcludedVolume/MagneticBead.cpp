#include "MagneticBead.h"

MagneticBead::MagneticBead(Chain * ch, arma::mat * bp_pos_backup, arma::cube * triads_backup, double bead_radius, double EV_radius, int attachment_bp ,const arma::vec & vec_triad2bead ) :
chain(ch),
BPS(*ch->get_BPS()),
bp_pos(ch->get_bp_pos()),
triads(ch->get_triads()),
bp_pos_backup(bp_pos_backup),
triads_backup(triads_backup),
num_bp(ch->get_num_bp()),
num_bps(ch->get_num_bps()),
disc_len(ch->get_disc_len()),
closed_topology(chain->topology_closed()),
bead_rad(bead_radius),
EV_rad(EV_radius),
attach_bp(attachment_bp),
counter(0),
counter_reject(0)
{
    min_dist = bead_rad+EV_rad;
    rel_vec  = triads->slice(attach_bp).t()*vec_triad2bead;
    rel_vec  = rel_vec / arma::norm(rel_vec) * (min_dist);
}


bool MagneticBead::check_bead_EV_double() {


    return false;
}

bool MagneticBead::check_bead_EV_current() {

    return false;
}



arma::vec MagneticBead::get_bead_pos() {
    return bp_pos->col(attach_bp) + triads->slice(attach_bp)*rel_vec;
}

arma::vec MagneticBead::get_bead_pos_backup() {
    return bp_pos_backup->col(attach_bp) + triads_backup->slice(attach_bp)*rel_vec;
}


double MagneticBead::doubleMove(const arma::vec & beadpos_past, const arma::vec & beadpos_curr, int bp_id) {
    arma::colvec p2,p2p,Delta_p,Delta_v,dvec;
    double lambda,dist,dist_primes;

    p2  = bp_pos_backup->col(bp_id);
    p2p = bp_pos       ->col(bp_id);

    dist_primes = arma::norm(p2p-beadpos_curr);
    if (dist_primes<min_dist) {
        return dist_primes;
    }
    Delta_p = beadpos_past-p2;
    Delta_v = beadpos_curr-p2p-Delta_p;

    // check if both were translated by the same vector
    dist = arma::dot(Delta_v,Delta_v);
    if (dist < 1e-10) {
        return dist_primes;
    }
    lambda = -arma::dot(Delta_p,Delta_v)/dist;
    // check if the closest approach is outside the moved interval
    if (lambda < 0) {
        return arma::norm(Delta_p);
    }
    if (lambda > 1) {
        return dist_primes;
    }
    dvec = Delta_p+lambda*Delta_v;
    dist = arma::norm(dvec);
    return dist;
}

