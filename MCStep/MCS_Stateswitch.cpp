#include "MCS_Stateswitch.h"

MCS_Stateswitch::MCS_Stateswitch(Chain * ch,const std::vector<long long int> & seedseq, int selrange_min, int selrange_max)
: MCStep(ch,seedseq)
{
    move_name = "MCS_Stateswitch";
    closed = chain->topology_closed();
    init_interface_energies();

    // set the range of the random generators selecting the hinges
    if (closed) {
        if (selrange_max >= num_bps/2)    selrange_max = num_bps/2;
    }
    // if (selrange_min < 2)            selrange_min = 2;
    if (selrange_max < selrange_min) selrange_max = selrange_min;

    if (selrange_min >= num_bps) {
        selrange_min = num_bps - 1;
    }
    if (selrange_max >= num_bps) {
        selrange_max = num_bps - 1;
    }

    if (closed) {
        decltype(genhinge1.param()) new_range(0, num_bps-1);
        genhinge1.param(new_range);

        rge_id1  = 0;
        rge_id2  = num_bps-1;
        rge_span = rge_id2-rge_id1;
    }
    else {
        int rfrom, rto;
        rfrom = 0;
        // if (!chain->fixed_first_orientation())  rfrom = 0;
        // else                                    rfrom = 1;
        rto   = num_bps-1-selrange_min;
        decltype(genhinge1.param()) new_range(rfrom, rto);
        genhinge1.param(new_range);

        rge_id1  = rfrom;
        rge_id2  = num_bps-1;
        rge_span = rge_id2-rge_id1;
    }
    decltype(genhingedist.param()) new_range(selrange_min, selrange_max);
    genhingedist.param(new_range);

////  OLD  //////////////////////////////////////////////////////////////////////
//    if (selrange_max >= num_bp/2)    selrange_max = num_bp/2;
//    if (selrange_min < 2)            selrange_min = 2;
//    if (selrange_max < selrange_min) selrange_max = selrange_min;
//
//    if (closed) {
//        decltype(genhinge1.param()) new_range(0, num_bp-1);
//        genhinge1.param(new_range);
//    }
//    else {
//        int rfrom, rto;
//        if (!chain->fixed_first_orientation())  rfrom = 0;
//        else                                    rfrom = 1;
//        rto   = num_bp-1-selrange_min;
//
//        decltype(genhinge1.param()) new_range(rfrom, rto);
//        genhinge1.param(new_range);
//    }
//    decltype(genhingedist.param()) new_range(selrange_min, selrange_max);
//    genhingedist.param(new_range);
////////////////////////////////////////////////////////////////////////////////

    // intervals of moved segments for exluded volume interactions
    // moved_intervals.push_back({-1,-1,1});
    // moved_intervals.push_back({0,-1,0});
    // moved_intervals.push_back({-1,num_bp-1,0});

    suitablility_key[MCS_FORCE_ACTIVE]              = true;  // force_active
    suitablility_key[MCS_TORQUE_ACTIVE]             = true;  // torque_active
    suitablility_key[MCS_TOPOLOGY_CLOSED]           = true;  // topology_closed
    suitablility_key[MCS_FIXED_LINK]                = true;  // link_fixed
    suitablility_key[MCS_FIXED_TERMINI]             = true;  // fixed_termini
    suitablility_key[MCS_FIXED_TERMINI_RADIAL]      = true;  // fixed_termini_radial
    suitablility_key[MCS_FIXED_FIRST_ORIENTATION]   = true;  // fixed_first_orientation
    suitablility_key[MCS_FIXED_LAST_ORIENTATION]    = true;  // fixed_last_orientation

    requires_constraint_check = false;
    requires_EV_check   = false;
    requires_pair_check = false;
}


MCS_Stateswitch::MCS_Stateswitch(Chain * ch,const std::vector<long long int> & seedseq, int selrange_min, int selrange_max, int range_id1, int range_id2)
: MCStep(ch,seedseq)
/*

    If the chain is open range_id2 cannot be larger than num_bp-2
*/
{
    move_name = "MCS_Stateswitch";
    closed = chain->topology_closed();
    init_interface_energies();

    // set the range of the random generators selecting the hinges
    if (closed) {
        if (selrange_max >= num_bps/2)    selrange_max = num_bps/2;
    }
    // if (selrange_min < 2)            selrange_min = 2;
    if (selrange_max < selrange_min) selrange_max = selrange_min;

    if (selrange_min >= num_bps) {
        selrange_min = num_bps - 1;
    }
    if (selrange_max >= num_bps) {
        selrange_max = num_bps - 1;
    }

    rge_id1  = pmod(range_id1  ,num_bps);
    rge_id2  = pmod(range_id2+1,num_bps);

    if (closed) {
        int range = pmod(rge_id2-rge_id1,num_bps);
        if (range == num_bps || range == 0) {
            rge_id1  = 0;
            rge_id2  = num_bps-1;
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
        if (rge_id1 < 0 || rge_id1 > num_bp-2) {
            std::cout << "Error invalid Range selection in MCS_Stateswitch. For open chains range_id1 needs to be within the range [0,num_bp-2]!" << std::endl;
            std::exit(0);
        }
        if (rge_id2 < 0 || rge_id2 > num_bp-1 || range_id2+1 == num_bp) {
            std::cout << "Error invalid Range selection in MCS_Stateswitch. For open chains range_id2 needs to be within the range [0,num_bp-2]!" << std::endl;
            std::exit(0);
        }
        if (rge_id2 < rge_id1) {
            std::cout << "Error invalid Range selection in MCS_Stateswitch. For open chains range_id2 needs to be larger or equal to range_id1!" << std::endl;
            std::exit(0);
        }
        if (chain->fixed_first_orientation() && rge_id1 == 0) {
            rge_id1 = 1;
        }

        int upper = larger(rge_id1,rge_id2-selrange_min);
        decltype(genhinge1.param()) new_range(rge_id1, upper);
        genhinge1.param(new_range);
    }
    rge_span = pmod(rge_id2-rge_id1,num_bps);

    decltype(genhingedist.param()) new_range(selrange_min, selrange_max);
    genhingedist.param(new_range);

    // intervals of moved segments for exluded volume interactions
    // moved_intervals.push_back({-1,-1,1});
    // moved_intervals.push_back({0,-1,0});
    // moved_intervals.push_back({-1,num_bp-1,0});

    suitablility_key[MCS_FORCE_ACTIVE]              = true;  // force_active
    suitablility_key[MCS_TORQUE_ACTIVE]             = true;  // torque_active
    suitablility_key[MCS_TOPOLOGY_CLOSED]           = true;  // topology_closed
    suitablility_key[MCS_FIXED_LINK]                = true;  // link_fixed
    suitablility_key[MCS_FIXED_TERMINI]             = true;  // fixed_termini
    suitablility_key[MCS_FIXED_TERMINI_RADIAL]      = true;  // fixed_termini_radial
    suitablility_key[MCS_FIXED_FIRST_ORIENTATION]   = true;  // fixed_first_orientation
    suitablility_key[MCS_FIXED_LAST_ORIENTATION]    = true;  // fixed_last_orientation

    requires_constraint_check = false;
    requires_EV_check   = true;
    requires_pair_check = true;
}



MCS_Stateswitch::~MCS_Stateswitch() {

}

void MCS_Stateswitch::update_settings() {
    // sigma = fac*sqrt(chain->get_T()/300)*std::sqrt(disc_len/0.34);
    init_interface_energies();
}

void MCS_Stateswitch::init_interface_energies() {
    std::cout << "Init interface_energies" << std::endl;
    interface_energies = arma::zeros(num_bps);
    for (unsigned i=0;i<BPS.size();i++) {
        interface_energies.at(i) = BPS[i]->bubble_interface_energy();
    }
}

bool MCS_Stateswitch::MC_move() {

    int             idA,idB;
    int             idAm1,idBp1;
    int             id;
    int             hingedist;
    double          theta;
    arma::mat       Rlab,TA_rot,TBn_rot;
	arma::colvec    Theta;
    double          deltaE = 0;
    int             new_state;

    idA       = genhinge1(gen);
    hingedist = genhingedist(gen);


    // std::cout << "idA       = " << idA << std::endl;
    // std::cout << "hingedist = " << hingedist << std::endl;

    idB       = idA + hingedist;
    if (closed) {
        idA    = pmod(idA,num_bps);
        idB    = pmod(idB,num_bps);
        idAm1  = pmod(idA-1,num_bps);
        idBp1  = pmod(idB+1,num_bps);

        if (closed_topol_restricted_range){
            int span = pmod(idB-rge_id1,num_bps);
            if (span > rge_span || span < hingedist) {
                // Check to keep idB within the selected range.
                idB = rge_id2;
                hingedist = pmod(idB-idA,num_bps);
            }
        }

    }
    else {
        if (idB > rge_id2 ) {
            idB = rge_id2;
            hingedist = idB-idA;
        }
        idAm1 = idA-1;
        idBp1 = idB+1;
    }

    for (int i=idA;i<=idA+hingedist;i++) {
        id = pmod(i,num_bps);
        if (states->at(id) == 1) {
            new_state = 2;
        }
        else {
            new_state = 1;
        }
        BPS[id]->propose_stateswitch(new_state);
    }

    for (int i=idA;i<=idA+hingedist;i++) {
        id = pmod(i,num_bps);
        deltaE += BPS[id]->eval_delta_energy();
    }

    ////////////////////////////////////
    // interface energy

    // std::cout << (*chain->get_states()).t() << std::endl;
    // std::cout << idAm1 << " " << idA << std::endl;

    if (idAm1 >= 0) {
        if (states->at(idAm1) == states->at(idA)) {
            deltaE += interface_energies(idAm1);
            // std::cout << interface_energies(idAm1) << std::endl;
        }
        else {
            deltaE -= interface_energies(idAm1);
            // std::cout << -interface_energies(idAm1) << std::endl;
        }
    }

    // std::cout << idB << " " << idBp1 << std::endl;
    if (idBp1 < num_bps) {
        if (states->at(idB) == states->at(idBp1)) {
            deltaE += interface_energies(idB);
            // std::cout << interface_energies(idAm1) << std::endl;
        }
        else {
            deltaE -= interface_energies(idB);
            // std::cout << -interface_energies(idAm1) << std::endl;
        }
    }

    ////////////////////////////////////

    // std::cout << "deltaE    = " << deltaE << std::endl;

    if (exp(-deltaE) <= uniformdist(gen)) {
        for (int i=idA;i<=idA+hingedist;i++) {
            id = pmod(i,num_bps);
            BPS[id]->set_switch(false);
        }
        // std::cout << (*chain->get_states()).t() << std::endl;
        // std::cout << arma::sum((*chain->get_states()))- chain->get_states()->size()) << std::endl;
		return false;
    }
    for (int i=idA;i<=idA+hingedist;i++) {
        id = pmod(i,num_bps);
        BPS[id]->set_switch(true);
        if (states->at(id) == 1) {
            new_state = 2;
        }
        else {
            new_state = 1;
        }
        states->at(id) = new_state;
    }
    // std::cout << (*chain->get_states()).t() << std::endl;
    // std::cout << arma::sum((*chain->get_states())) - chain->get_states()->size() << std::endl;
    return true;
}




