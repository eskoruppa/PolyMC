#include "MCS_ClusterTwist.h"


/*
    For open chains this move will never touch the termini
*/

MCS_ClusterTwist::MCS_ClusterTwist(Chain * ch,const std::vector<long long int> & seedseq, int selrange_min, int selrange_max)
: MCStep(ch,seedseq)
{
    move_name = "MCS_ClusterTw";
    /*
    Initializes the controling variables of the ClusterTwist move
    */

    /*
        fac for plasmids
    */
    sigfac = 3;
    arma::mat covmat = chain->get_avg_cov();
    sigma = sigfac*sqrt(covmat(2,2));

    closed = chain->topology_closed();

    // set the range of the random generators selecting the hinges
    if (selrange_max >= num_bp/2)    selrange_max = num_bp/2;
    if (selrange_min < 2)            selrange_min = 2;
    if (selrange_max < selrange_min) selrange_max = selrange_min;

    if (closed) {
        decltype(genhinge1.param()) new_range(0, num_bp-1);
        genhinge1.param(new_range);

        rge_id1  = 0;
        rge_id2  = num_bp-1;
        rge_span = rge_id2-rge_id1;
    }
    else {
        decltype(genhinge1.param()) new_range(1, num_bp-2-selrange_min);
        genhinge1.param(new_range);

        rge_id1  = 1;
        rge_id2  = num_bp-2;
        rge_span = rge_id2-rge_id1;
    }
    decltype(genhingedist.param()) new_range(selrange_min, selrange_max);
    genhingedist.param(new_range);

//    // intervals of moved segments for exluded volume interactions
//    moved_intervals.push_back({-1,-1,1});
//    moved_intervals.push_back({0,-1,0});
//    moved_intervals.push_back({-1,num_bp-1,0});

    /*
        Quick evaluation is posible for local isotropic hamiltonians. Requires elastic energy evaluation
        only at the hinge site and allows for triad rotations to be carried out only if move is accepted.
    */
    if (chain->local_isotropic_hamiltonian()) {
        quick_eval = true;
    }
    else {
        quick_eval = false;
    }

    suitablility_key[MCS_FORCE_ACTIVE]              = true;  // force_active
    suitablility_key[MCS_TORQUE_ACTIVE]             = true;  // torque_active
    suitablility_key[MCS_TOPOLOGY_CLOSED]           = true;  // topology_closed
    suitablility_key[MCS_FIXED_LINK]                = true;  // link_fixed
    suitablility_key[MCS_FIXED_TERMINI]             = true;  // fixed_termini
    suitablility_key[MCS_FIXED_TERMINI_RADIAL]      = true;  // fixed_termini_radial
    suitablility_key[MCS_FIXED_FIRST_ORIENTATION]   = true;  // fixed_first_orientation
    suitablility_key[MCS_FIXED_LAST_ORIENTATION]    = true;  // fixed_last_orientation

    requires_EV_check=false;
}

MCS_ClusterTwist::MCS_ClusterTwist(Chain * ch,const std::vector<long long int> & seedseq, int selrange_min, int selrange_max, int range_id1, int range_id2)
: MCStep(ch,seedseq)
{
    move_name = "MCS_ClusterTw";
    /*
    Initializes the controling variables of the ClusterTwist move
    */

    /*
        fac for plasmids
    */
    sigfac = 3;
    arma::mat covmat = chain->get_avg_cov();
    sigma = sigfac*sqrt(covmat(2,2));

    closed = chain->topology_closed();

    // set the range of the random generators selecting the hinges
    if (selrange_max >= num_bp/2)    selrange_max = num_bp/2;
    if (selrange_min < 2)            selrange_min = 2;
    if (selrange_max < selrange_min) selrange_max = selrange_min;

    rge_id1  = pmod(range_id1  ,num_bp);
    rge_id2  = pmod(range_id2+1,num_bp);

    if (closed) {
        int range = pmod(rge_id2-rge_id1,num_bp);
        if (range == num_bp || range == 0) {
            rge_id1  = 0;
            rge_id2  = num_bp-1;
            decltype(genhinge1.param()) new_range(rge_id1,rge_id2);
            genhinge1.param(new_range);
        }
        else {
            closed_topol_restricted_range = true;
            int upper = larger(rge_id1,rge_id1+range-selrange_min);
            decltype(genhinge1.param()) new_range(rge_id1, upper);
            genhinge1.param(new_range);
        }
    }
    else {
        if (rge_id1 < 1 || rge_id1 > num_bp-2) {
            std::cout << "Error invalid Range selection in MCS_ClusterTwist. For open chains range_id1 needs to be within the range [1,num_bp-2]!" << std::endl;
            std::exit(0);
        }
        if (rge_id2 < 1 || rge_id2 > num_bp-1) {
            std::cout << "Error invalid Range selection in MCS_ClusterTwist. For open chains range_id2 needs to be within the range [1,num_bp-2]!" << std::endl;
            std::exit(0);
        }
        if (rge_id2 < rge_id1) {
            std::cout << "Error invalid Range selection in MCS_ClusterTwist. For open chains range_id2 needs to be larger or equal to range_id1!" << std::endl;
            std::exit(0);
        }

        int upper = larger(rge_id1,rge_id2-selrange_min);
        decltype(genhinge1.param()) new_range(rge_id1, upper);
        genhinge1.param(new_range);
    }
    rge_span = pmod(rge_id2-rge_id1,num_bp);

    decltype(genhingedist.param()) new_range(selrange_min, selrange_max);
    genhingedist.param(new_range);

//    // intervals of moved segments for exluded volume interactions
//    moved_intervals.push_back({-1,-1,1});
//    moved_intervals.push_back({0,-1,0});
//    moved_intervals.push_back({-1,num_bp-1,0});

    /*
        Quick evaluation is posible for local isotropic hamiltonians. Requires elastic energy evaluation
        only at the hinge site and allows for triad rotations to be carried out only if move is accepted.
    */
    if (chain->local_isotropic_hamiltonian()) {
        quick_eval = true;
    }
    else {
        quick_eval = false;
    }

    suitablility_key[MCS_FORCE_ACTIVE]              = true;  // force_active
    suitablility_key[MCS_TORQUE_ACTIVE]             = true;  // torque_active
    suitablility_key[MCS_TOPOLOGY_CLOSED]           = true;  // topology_closed
    suitablility_key[MCS_FIXED_LINK]                = true;  // link_fixed
    suitablility_key[MCS_FIXED_TERMINI]             = true;  // fixed_termini
    suitablility_key[MCS_FIXED_TERMINI_RADIAL]      = true;  // fixed_termini_radial
    suitablility_key[MCS_FIXED_FIRST_ORIENTATION]   = true;  // fixed_first_orientation
    suitablility_key[MCS_FIXED_LAST_ORIENTATION]    = true;  // fixed_last_orientation

    requires_EV_check=false;
}


MCS_ClusterTwist::~MCS_ClusterTwist() {

}

void MCS_ClusterTwist::update_settings() {
    arma::mat covmat = chain->get_avg_cov();
    sigma = sigfac*sqrt(covmat(2,2));
}

bool MCS_ClusterTwist::MC_move() {

    int         idA,idB;
    int         idAm1,idBm1;
    int         hingedist;
    double      deltaE = 0;
    double      dtheta;
    arma::mat   Rz;

    double choice = uniformdist(gen);

    idA       = genhinge1(gen);
    hingedist = genhingedist(gen);
    idB       = idA + hingedist;
    if (closed) {
        idA    = pmod(idA,num_bp);
        idAm1  = pmod(idA-1,num_bp);
        idB    = pmod(idB,num_bp);
        if (closed_topol_restricted_range){
            int span = pmod(idB-rge_id1,num_bp);
            if (span > rge_span || span < hingedist) {
                // Check to keep idBn within the selected range.
                idB = rge_id2;
                hingedist = pmod(idB-idA,num_bp);
            }
        }
        idBm1  = pmod(idB-1,num_bp);
    }
    else {
        idAm1 = idA-1;
        if (idB > rge_id2 ) {
            idB = rge_id2;
            hingedist = idB-idA;
        }
        idBm1 = idB-1;
    }

    dtheta = normaldist(gen)*sigma;
    Rz     = Rotz(dtheta);

    if (quick_eval) {
        /*
            If the hamiltonian is local and bending is isotropic the energy will only be evaluated at the hinge
            and at the terminal bead if the torsional trap or an external torque constraint is active.
        */

        arma::mat TA, TBm1;

        TA   = triads->slice(idA  )*Rz;
        TBm1 = triads->slice(idBm1)*Rz;

        BPS[idAm1]->propose_move(triads->slice(idAm1), TA);
        BPS[idBm1]->propose_move(TBm1, triads->slice(idB));

        deltaE += BPS[idAm1]->eval_delta_energy();
        deltaE += BPS[idBm1]->eval_delta_energy();

        changed_bps = {idAm1,idBm1};

        if (exp(-deltaE) <= choice) {
            return false;
        }

        triads->slice(idA)   = TA;
        triads->slice(idBm1) = TBm1;
        int id;
        for (int i=idA+1;i<idA+hingedist-1;i++) {
            id = pmod(i,num_bp);
            triads->slice(id) = triads->slice(id)*Rz;
        }
        return true;
    }
    else {

        arma::cube triads_cp = *triads;
        changed_bps = arma::zeros<arma::ivec>(hingedist+1);

        int id, idm1;
        id = pmod(idA-1,num_bp);
        for (int i=0;i<hingedist;i++) {
            idm1 = id;
            id   = pmod(idm1+1,num_bp);
            triads_cp.slice(id) = triads_cp.slice(id)*Rz;
            BPS[idm1]->propose_move(triads_cp.slice(idm1), triads_cp.slice(id));
            changed_bps(i) = idm1;
        }
        BPS[idBm1]->propose_move(triads_cp.slice(idBm1), triads_cp.slice(idB));
        changed_bps(hingedist) = idBm1;

        for (int i=idAm1;i<=idAm1+hingedist;i++) {
            id = pmod(i,num_bp);
            deltaE += BPS[id]->eval_delta_energy();
        }

        if (exp(-deltaE) <= choice) {
            return false;
        }

        *triads = triads_cp;
        return true;
    }
    return true;
}
