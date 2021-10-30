#include "MCS_Reptation_Thermal.h"


MCS_RepTherm::MCS_RepTherm(Chain * ch, int rep_size_min, int rep_size_max, int rep_dist_min, int rep_dist_max, bool twist_translation_active,bool thermalize)
: MCStep(ch),
rep_size_min(rep_size_min),
rep_size_max(rep_size_max),
rep_dist_min(rep_dist_min),
rep_dist_max(rep_dist_max),
twist_translation_active(twist_translation_active),
thermalize(thermalize)
{
    move_name="MCS_RepTherm";

    MC_thermal_rot     = new MCS_CSrot(chain,rep_size_min,rep_size_max);
    MC_pre_thermal_rot = new MCS_CSrot(chain,rep_size_min,rep_size_max);

    closed = chain->topology_closed();

    /*
        Repations should not span over more than half the length of the chain
    */
    if (rep_dist_max>num_bp/2) {
        rep_dist_max=num_bp/2;
        if (rep_dist_min>rep_dist_max) {
            rep_dist_min = rep_dist_max;
        }
    }

    if (2*(LEVER_SIZE_FAC+1)*rep_size_max>=num_bp/2) {
        rep_size_max = num_bp/(4*(LEVER_SIZE_FAC+1))-1;
        if (rep_size_min>rep_size_max) {
            rep_size_min=rep_size_max;
        }
    }


    if (closed) {
        decltype(gen_idFirst.param()) new_range_idFirst(0, num_bp-1);
        gen_idFirst.param(new_range_idFirst);

    }
    else {
        int min_idFirst = THERMALIZE_DIST+1;
        /*
            The +2 is added so that the first triad does not have to be rotated
        */
        int max_idFirst = num_bp-THERMALIZE_DIST-2 - rep_dist_min - (2*LEVER_SIZE_FAC+1)*rep_size_min;

        decltype(gen_idFirst.param()) new_range_idFirst(min_idFirst, max_idFirst);
        gen_idFirst.param(new_range_idFirst);
    }

    decltype(gen_repsize.param()) new_range_repsize(rep_size_min, rep_size_max);
    gen_repsize.param(new_range_repsize);

    decltype(gen_repdist.param()) new_range_repdist(rep_dist_min, rep_dist_max);
    gen_repdist.param(new_range_repdist);

//    // intervals of moved segments for exluded volume interactions
    moved_intervals.push_back({-1,-1,2000});
    moved_intervals.push_back({-1,-1,2000});
    moved_intervals.push_back({-1,-1,2000});
    moved_intervals.push_back({-1,-1,2000});
    moved_intervals.push_back({-1,-1,2000});
    /*
    Typically there will be 2 moved intervals. However one can be split in
    2 intervals when it crosses the boundary on a closed chain
    */
    suitablility_key[MCS_FORCE_ACTIVE]              = true;  // force_active
    suitablility_key[MCS_TORQUE_ACTIVE]             = true;  // torque_active
    suitablility_key[MCS_TOPOLOGY_CLOSED]           = true;  // topology_closed
    suitablility_key[MCS_FIXED_LINK]                = true;  // link_fixed
    suitablility_key[MCS_FIXED_TERMINI]             = true;  // fixed_termini
    suitablility_key[MCS_FIXED_TERMINI_RADIAL]      = true;  // fixed_termini_radial
    suitablility_key[MCS_FIXED_FIRST_ORIENTATION]   = true;  // fixed_first_orientation
    suitablility_key[MCS_FIXED_LAST_ORIENTATION]    = true;  // fixed_last_orientation

    requires_EV_check=true;
}

MCS_RepTherm::~MCS_RepTherm() {
    delete MC_thermal_rot;
}


void MCS_RepTherm::update_settings() {
    MC_thermal_rot->update_settings();
}

bool MCS_RepTherm::MC() {
    bool accepted = MC_move();
    if (accepted) {
        count_accept++;
    }
    return accepted;
}

void MCS_RepTherm::set_excluded_volume(ExVol* EVol) {
    ev_active = true;
    EV        = EVol;
    MC_pre_thermal_rot->set_excluded_volume(EV);
//    MC_thermal_rot->set_excluded_volume(EV);
}


bool MCS_RepTherm::MC_move() {

    double E_init, E_final;

    arma::colvec rem_v1,rem_v2, rem_w;
    double len_rem_v1,len_rem_v2,len_rem_w;

    /*
        Determine forward or backward reptation
    */

    if (randbinary(gen)==1) rep_forward=true;
    else                    rep_forward=false;

    /*
        Determine hinges
    */

    int trial_count=0;
    bool found = false;
    while (trial_count<CAP_RETRIES) {
        trial_count++;

        // generate reptation segment size and reptation distance
        rep_size    = gen_repsize(gen);
        rep_dist    = gen_repdist(gen);

        // size of a rotated segment
        lever_size  = LEVER_SIZE_FAC*rep_size;
        // size of the untouched segment in the middle
        intseg_size = rep_dist-2*lever_size;
        if (intseg_size<MIN_INTSEG_SIZE*lever_size) {
            intseg_size = MIN_INTSEG_SIZE*lever_size;
            rep_dist    = intseg_size + 2*lever_size;
        }
        // total size from the left hinge of the first segment to the right hinge
        // of the last segment
        total_size  = rep_dist + 2*lever_size + rep_size;

        // generate the left hinge of the first segment
        idFirst = gen_idFirst(gen);

        if (!closed) {
            int overshoot = idFirst+total_size - num_bp + THERMALIZE_DIST+2;
            if (overshoot > 0) {
                idFirst = idFirst-overshoot;
            }
        }

        /*
            Calculate all the hinges before and after the move
            TODO: For open chains the pmod can be omitted.
        */
        if (rep_forward) {
            rem_hA  = idFirst;
            rem_hv1 = pmod(rem_hA+lever_size , num_bp);
            rem_hv2 = pmod(rem_hv1+rep_size   , num_bp);
            rem_hB  = pmod(rem_hv2+lever_size , num_bp);

            rem_hAp = rem_hA;
            rem_hBp = pmod(rem_hB-rep_size, num_bp);
            rem_hvp = pmod(rem_hv1         , num_bp);

            add_hA  = pmod(rem_hB + intseg_size , num_bp);
            add_hB  = pmod(rem_hA + total_size  , num_bp);
            add_hv  = pmod(add_hA + lever_size  , num_bp);

            add_hAp  = pmod(add_hA-rep_size , num_bp);
            add_hBp  = add_hB;
            add_hv2p = add_hv;
            add_hv1p = pmod(add_hv-rep_size  , num_bp);

            idLast = add_hB;
        }
        else {
            add_hA  = idFirst;
            add_hv  = pmod(add_hA + lever_size , num_bp);
            add_hB  = pmod(add_hv  + lever_size , num_bp);

            add_hAp  = add_hA;
            add_hBp  = pmod(add_hB+rep_size , num_bp);
            add_hv1p = add_hv;
            add_hv2p = pmod(add_hv1p+rep_size  , num_bp);;

            rem_hA  = pmod(add_hB+intseg_size , num_bp);
            rem_hv1 = pmod(rem_hA+lever_size  , num_bp);
            rem_hv2 = pmod(rem_hv1+rep_size    , num_bp);
            rem_hB  = pmod(rem_hv2+lever_size  , num_bp);

            rem_hAp = pmod(rem_hA+rep_size, num_bp);
            rem_hBp = rem_hB;
            rem_hvp = rem_hv2;

            idLast = rem_hB;
        }




        /*
            Check if this move is possible
            Possibly repeat the generation of the hinges until a possible move is found (repetitions capped)
        */

        rem_w  = pos->col(rem_hB)-pos->col(rem_hA);
        len_rem_w  = arma::norm(rem_w);
//        rem_v1 = pos->col(rem_hv1)-pos->col(rem_hA);
//        rem_v2 = pos->col(rem_hB)-pos->col(rem_hv2);
//        len_rem_v1 = arma::norm(rem_v1);
//        len_rem_v2 = arma::norm(rem_v2);
//        if ((len_rem_v1+len_rem_v2)>len_rem_w) {
        if (SELECT_EXTRA_CRITERION_FACTOR*2*lever_size*disc_len>len_rem_w) {
            found=true;
            break;
        }
    }

    if (found==false) {
        return false;
    }

//    Pre_Thermalize();

//    cout << "#####################" << endl;
//    cout << endl;
//    cout << "rem_hA   " << rem_hA << endl;
//    cout << "rem_hv1  " << rem_hv1 << endl;
//    cout << "rem_hv2  " << rem_hv2 << endl;
//    cout << "rem_hB   " << rem_hB << endl;
//    cout << "add_hA   " << add_hA << endl;
//    cout << "add_hv   " << add_hv << endl;
//    cout << "add_hB   " << add_hB << endl;
//    cout << "---" << endl;
//    cout << "rem_hAp  " << rem_hAp << endl;
//    cout << "rem_hvp  " << rem_hvp << endl;
//    cout << "rem_hBp  " << rem_hBp << endl;
//    cout << "add_hAp  " << add_hAp << endl;
//    cout << "add_hv1p " << add_hv1p << endl;
//    cout << "add_hv2p " << add_hv2p << endl;
//    cout << "add_hBp  " << add_hBp << endl;
//    cout << endl;




    idEnergyFirst = pmod(idFirst-THERMALIZE_DIST-1,num_bp);
    idEnergyLast  = pmod(idLast +THERMALIZE_DIST-1,num_bp);
    idThermFirst  = pmod(idFirst-THERMALIZE_DIST  ,num_bp);
    idThermLast   = pmod(idLast +THERMALIZE_DIST  ,num_bp);

    int E_rem_from = pmod(rem_hAp-THERMALIZE_DIST-1,num_bp);
    int E_rem_to   = pmod(rem_hBp+THERMALIZE_DIST-1,num_bp);
    int E_add_from = pmod(add_hAp-THERMALIZE_DIST-1,num_bp);
    int E_add_to   = pmod(add_hBp+THERMALIZE_DIST-1,num_bp);

    double E_rem_init;
    double E_rem_final;
    double E_add_init;
    double E_add_final;


    /*
        Calculate energetic contribution before the Move
    */
//    chain->extract_energy_select(idEnergyFirst,idEnergyLast);
    E_init = chain->extract_energy(idEnergyFirst,idEnergyLast);


//    chain->extract_energy_select(E_rem_from,E_rem_to);
//    E_rem_init = chain->extract_energy(E_rem_from,E_rem_to);
//
//    chain->extract_energy_select(E_add_from,E_add_to);
//    E_add_init = chain->extract_energy(E_add_from,E_add_to);

    /*
        Make Backup of relevant configuration
    */

    conf_make_backup(idThermFirst,idThermLast);

    /*
        Calculate Twist translation
    */

    if (twist_translation_active) {
        translate_twist();
    }


    bool acpt;
    /*
        REM Move
    */
    acpt = REM_move();
    if (!acpt) {
        conf_revert_to_backup(idThermFirst,idThermLast);
        return false;
    }
    /*
        MID Move
    */
    MID_move();
    /*
        ADD Move
    */
    acpt = ADD_move();
    if (!acpt) {
        conf_revert_to_backup(idThermFirst,idThermLast);
        return false;
    }


    chain->cal_energy_propose(idEnergyFirst,idEnergyLast);
    chain->cal_energy_eval(idEnergyFirst,idEnergyLast);
    chain->cal_energy_set_move(idEnergyFirst,idEnergyLast,true);

    Thermalize_Hinges();
//    Thermalize();

    /*
        Calculate New Energy
    */
//    chain->cal_energy_propose(E_rem_from,E_rem_to);
//    E_rem_final = chain->cal_energy_eval(E_rem_from,E_rem_to);
//
//    chain->cal_energy_propose(E_add_from,E_add_to);
//    E_add_final = chain->cal_energy_eval(E_add_from,E_add_to);


    chain->cal_energy_propose(idEnergyFirst,idEnergyLast);

    E_final = chain->cal_energy_eval(idEnergyFirst,idEnergyLast);

    double deltaE = E_final - E_init;



    if (deltaE>1000 && false) {


//        cout << "E_rem_init  = " << E_rem_init  << endl;
//        cout << "E_rem_final = " << E_rem_final << endl;
//        cout << "E_add_init  = " << E_add_init  << endl;
//        cout << "E_add_final = " << E_add_final << endl;
//        cout << "E_init      = " << E_init << endl;
//        cout << "E_final     = " << E_final << endl;
//        cout << "delta E     = " << deltaE << endl;
//        cout << "<delta E>   = " << avg_deltaE/avg_deltaE_count << endl;
//        cout << endl;

        arma::mat M0 = *BPS[0]->get_M0();

        double check_E = 0;
        arma::colvec cE(3,arma::fill::zeros);

        arma::colvec Theta;
        int id1,id2,start;

        if (rep_forward) start = rem_hA-THERMALIZE_DIST;
        else             start = add_hA-THERMALIZE_DIST;
//        cout << "###########################" << endl;
//        cout << "###########################" << endl;
        for (int i =0;i<=total_size+2*THERMALIZE_DIST;i++) {
            id1 = pmod(start-1+i,num_bp);
            id2 = pmod(start+i,num_bp);
//            Theta = ExtractTheta(triads->slice(id1).t()*triads->slice(id2));
            Theta = ExtractTheta((*BPS[id1]->get_R0()).t()* triads->slice(id1).t()*triads->slice(id2));

            if (Theta.has_nan()) {
                cout << triads->slice(id1);
                cout << triads->slice(id2);
            }

//            if (id2==rem_hA) cout << "rem_hA" << endl;
//            if (id2==rem_hAp) cout << "rem_hAp" << endl;
//            if (id2==rem_hvp) cout << "rem_hvp" << endl;
//            if (id2==rem_hv1) cout << "rem_hv1" << endl;
//            if (id2==rem_hv2) cout << "rem_hv2" << endl;
//            if (id2==rem_hB) cout << "rem_hB" << endl;
//            if (id2==rem_hBp) cout << "rem_hBp" << endl;
//
//            if (id2==add_hA) cout << "add_hA" << endl;
//            if (id2==add_hAp) cout << "add_hAp" << endl;
//            if (id2==add_hv) cout << "add_hv" << endl;
//            if (id2==add_hv1p) cout << "add_hv1p" << endl;
//            if (id2==add_hv2p) cout << "add_hv2p" << endl;
//            if (id2==add_hB) cout << "add_hB" << endl;
//            if (id2==add_hBp) cout << "add_hBp" << endl;

//            cout << id2 << "  ---------------------" << endl;
//            cout << Theta.t()*180/M_PI ;

            double oneE;
            Theta = ExtractTheta((*BPS[id1]->get_R0()).t()* triads->slice(id1).t()*triads->slice(id2));
            for (int d=0;d<3;d++) {
                oneE = Theta(d)*Theta(d)*M0(d,d);
//                cout << oneE << " ";
                cE(d) += oneE;
            }
//            cout << endl;

            check_E += arma::dot(Theta, M0 * Theta);

        }
        cout << "##################################" << endl;
        cout << "<delta E>   = " << avg_deltaE/avg_deltaE_count << endl;
        cout << "delta E     = " << deltaE  << endl;
        cout << "delta E_rem = " << E_rem_final-E_rem_init  << endl;
        cout << "delta E_add = " << E_add_final-E_add_init  << endl;
        cout << "check new E = " << E_final-check_E  << endl;

        cout << "-----------------" << endl;
        cout << "New     energy = " << check_E << endl;
        cout << "Initial energy = " << E_init << endl;
        cout << cE.t();
    }

    if (deltaE!=deltaE) {
        cout << "Energy is NaN!" << endl;
        std::exit(0);
        return false;
    }
    else {
        avg_deltaE      += deltaE;
        avg_deltaE_count++;
    }

//    cout << "<delta E>   = " << avg_deltaE/avg_deltaE_count << endl;


    /*
         Metropolis Step
    */


//    if  (exp(-deltaE) <= uniformdist(gen)){
//        conf_revert_to_backup(idThermFirst,idThermLast);
//		return false;
//    }

    deltaE = deltaE/std::sqrt((6*THERMALIZE_DIST)+rep_size-6);
    deltaE = deltaE/((6*THERMALIZE_DIST)+rep_size-6);

    count_step++;
    bool accepted (exp(-deltaE) > uniformdist(gen));
    if (accepted) {
//        cout << "deltaE   = " << deltaE << endl;
        if (ev_active) {
            accepted = Check_EV();
//            if (!accepted) {
//                cout << "EV violation!" << endl;
//                accepted = true;
//            }
        }
        chain->cal_energy_set_move(idEnergyFirst,idEnergyLast,accepted);
    }
    else {
        conf_revert_to_backup(idThermFirst,idThermLast);
        chain->cal_energy_set_move(idEnergyFirst,idEnergyLast,false);
    }

    if (!accepted) {
        chain->cal_energy_propose(idEnergyFirst,idEnergyLast);
        chain->cal_energy_eval(idEnergyFirst,idEnergyLast);
        chain->cal_energy_set_move(idEnergyFirst,idEnergyLast,true);
    }


//    conf_revert_to_backup(idThermFirst,idThermLast);
//    if (accepted) {
//        cout << "------------------------------" << endl;
//        cout << "rep size:        " << rep_size << endl;
//        cout << "<delta E>      = " << avg_deltaE/avg_deltaE_count << endl;
//        cout << "delta E        = " << deltaE  << endl;
//        cout << "New     energy = " << E_final << endl;
//        cout << "Initial energy = " << E_init << endl;
//
//        conf_revert_to_backup(idThermFirst,idThermLast);
//
//        int N = 10000;
//
//        auto timer_start = std::chrono::high_resolution_clock::now();
//        auto timer_finish = std::chrono::high_resolution_clock::now();
//        chrono::duration<double> timer_elapsed;
//
//        for (int L=0;L<N;L++) {
//            REM_move();
//            conf_revert_to_backup(idThermFirst,idThermLast);
//        }
//
//        timer_finish = std::chrono::high_resolution_clock::now();
//        timer_elapsed       = timer_finish - timer_start;
//        cout << "elapsed time REM:   " << timer_elapsed.count() << " s\n";
//        timer_start = std::chrono::high_resolution_clock::now();
//
//        REM_move();
//        MID_move();
//
//        for (int L=0;L<N;L++) {
//            ADD_move();
//            conf_revert_to_backup(idThermFirst,idThermLast);
//        }
//
//        timer_finish = std::chrono::high_resolution_clock::now();
//        timer_elapsed       = timer_finish - timer_start;
//        cout << "elapsed time ADD:   " << timer_elapsed.count() << " s\n";
//        timer_start = std::chrono::high_resolution_clock::now();
//
//        for (int L=0;L<N;L++) {
//            chain->cal_energy_propose(idEnergyFirst,idEnergyLast);
//            chain->cal_energy_eval(idEnergyFirst,idEnergyLast);
//            chain->cal_energy_set_move(idEnergyFirst,idEnergyLast,true);
//        }
//
//        timer_finish = std::chrono::high_resolution_clock::now();
//        timer_elapsed       = timer_finish - timer_start;
//        cout << "elapsed time EN1:   " << timer_elapsed.count() << " s\n";
//        timer_start = std::chrono::high_resolution_clock::now();
//
//
//        for (int L=0;L<N;L++) {
//            arma::mat M0 = *BPS[0]->get_M0();
//            arma::mat R0T = (*BPS[0]->get_R0()).t();
//            double check_E = 0;
//            arma::colvec Theta;
//            int id1,id2,start;
//
//            if (rep_forward) start = rem_hA-THERMALIZE_DIST;
//            else             start = add_hA-THERMALIZE_DIST;
//            for (int i =0;i<=total_size+2*THERMALIZE_DIST;i++) {
//                id1 = pmod(start-1+i,num_bp);
//                id2 = pmod(start+i,num_bp);
//
//                Theta = ExtractTheta(R0T* triads->slice(id1).t()*triads->slice(id2));
//                check_E += arma::dot(Theta, M0 * Theta);
//            }
//        }
//
//        timer_finish = std::chrono::high_resolution_clock::now();
//        timer_elapsed       = timer_finish - timer_start;
//        cout << "elapsed time EN2:   " << timer_elapsed.count() << " s\n";
//        timer_start = std::chrono::high_resolution_clock::now();
//
//        for (int L=0;L<N/10;L++) {
//            Thermalize();
//            conf_revert_to_backup(idThermFirst,idThermLast);
//        }
//        timer_finish = std::chrono::high_resolution_clock::now();
//        timer_elapsed       = timer_finish - timer_start;
//        cout << "elapsed time THA:   " << timer_elapsed.count() << " s\n";
//        timer_start = std::chrono::high_resolution_clock::now();
//
//        for (int L=0;L<N;L++) {
//            Thermalize_Hinges();
//            conf_revert_to_backup(idThermFirst,idThermLast);
//        }
//        timer_finish = std::chrono::high_resolution_clock::now();
//        timer_elapsed       = timer_finish - timer_start;
//        cout << "elapsed time THH:   " << timer_elapsed.count() << " s\n";
//
//    }


    return accepted;
}


void MCS_RepTherm::conf_make_backup(int idThermFirst, int idThermLast) {

//    if (idThermFirst>idThermLast) {
//        backup_bp_pos1 = pos->cols(idThermFirst,num_bp-1);
//        backup_triads1 = triads->slices(idThermFirst,num_bp-1);
//
//        backup_bp_pos2 = pos->cols(0,idThermLast);
//        backup_triads2 = triads->slices(0,idThermLast);
//    }
//    else {
//        backup_bp_pos1 = pos->cols(idThermFirst,idThermLast);
//        backup_triads1 = triads->slices(idThermFirst,idThermLast);
//    }
    backup_bp_pos1 = *pos;
    backup_triads1 = *triads;
}

void MCS_RepTherm::conf_revert_to_backup(int idThermFirst, int idThermLast) {

//    if (idThermFirst>idThermLast) {
//        pos->cols(idThermFirst,num_bp-1)      = backup_bp_pos1;
//        triads->slices(idThermFirst,num_bp-1) = backup_triads1;
//
//        pos->cols(0,idThermLast)      = backup_bp_pos2;
//        triads->slices(0,idThermLast) = backup_triads2;
//    }
//    else {
//        pos->cols(idThermFirst,idThermLast)      = backup_bp_pos1;
//        triads->slices(idThermFirst,idThermLast) = backup_triads1;
//    }
    *pos    = backup_bp_pos1;
    *triads = backup_triads1;
    EV->set_backup_conf(&backup_bp_pos1,&backup_triads1);
}


void MCS_RepTherm::translate_twist() {
    double twist_translation=0;
    int id;
    arma::mat Rtw;
    for (int i=0;i<rep_size;i++) {
        id = pmod(rem_hv1+i,num_bp);
        twist_translation += BPS[id]->get_T0_component(2);
    }
    if (rep_forward) {
        Rtw = Rotz(-twist_translation);
        for (int i=0;i<rep_dist;i++) {
            id = pmod(rem_hv2+i,num_bp);
            triads->slice(id) = triads->slice(id)*Rtw;
        }
    }
    else {
        Rtw = Rotz(twist_translation);
        for (int i=0;i<rep_dist;i++) {
            id = pmod(add_hv+i,num_bp);
            triads->slice(id) = triads->slice(id)*Rtw;
        }
    }
}



bool MCS_RepTherm::REM_move() {
/*
    TODO remove double entry member variables (e.g. rem_v1 ..)
*/

    arma::colvec rem_v1,rem_v2,rem_v3, rem_w;
    double       len_rem_v1,len_rem_v2,len_rem_w;
    double       len_v1_orth, len_v2_orth;
    arma::colvec x,y,z;
    arma::colvec nv2_orth;
    double theta,phi;

//    cout << "##########" << endl;
//    cout << (pos->col(rem_hv2)-pos->col(rem_hv1)).t();

    arma::mat rot_align1, rot_align2, rot_hightalign;

    rem_v1      = pos->col(rem_hv1)-pos->col(rem_hA);
    rem_v2      = pos->col(rem_hB)-pos->col(rem_hv2);
    rem_v3      = pos->col(rem_hv2)-pos->col(rem_hA);
    rem_w       = pos->col(rem_hB) -pos->col(rem_hA);
    len_rem_v1  = arma::norm(rem_v1);
    len_rem_v2  = arma::norm(rem_v2);
    len_rem_w   = arma::norm(rem_w);

    x           = rem_w/len_rem_w;
    y           = rem_v1-arma::dot(rem_v1,x)*x;
    len_v1_orth = arma::norm(y);
    y           = y/len_v1_orth;
    z           = arma::cross(x,y);
    nv2_orth    = rem_v3-arma::dot(rem_v3,x)*x;
    len_v2_orth = arma::norm(nv2_orth);
    nv2_orth    = nv2_orth/len_v2_orth;

    // align the orthogoal components of v1 and v2
    theta = std::atan2(arma::dot(nv2_orth,z),arma::dot(nv2_orth,y));
    rot_align1 = getRotMat(0.5*theta*x);
    rot_align2 = getRotMat(-0.5*theta*x);

    // rotate y and z
    y = rot_align1*y;
    z = rot_align1*z;

    // equalize the orthogonal componnts of v1 and v2
    phi   = std::atan(    (len_v2_orth-len_v1_orth )
                        / ( std::sqrt(len_rem_v1*len_rem_v1-len_v1_orth*len_v1_orth)
                           +std::sqrt(len_rem_v2*len_rem_v2-len_v2_orth*len_v2_orth)));
    rot_hightalign = getRotMat(phi*z);

    // calculate total rotation
    rot_align1 = rot_hightalign*rot_align1;
    rot_align2 = rot_hightalign*rot_align2;


    /*
        Calculate Orthogonal and Parallel components of the tangent vectors
    */
    int id_from, id_to;
    arma::mat seg1_para(3,lever_size);
    arma::mat seg1_orth(3,lever_size);
    arma::mat seg2_para(3,lever_size);
    arma::mat seg2_orth(3,lever_size);

    arma::colvec tangent;
    arma::vec t_orthsq(2*lever_size);
    for (int i=1;i<=lever_size;i++) {
        id_from = pmod(rem_hv1-i,num_bp);
        id_to   = pmod(rem_hvp-i,num_bp);
        triads->slice(id_to) = rot_align1*triads->slice(id_from);

        tangent = triads->slice(id_to).col(2)*disc_len;
        seg1_para.col(lever_size-i) = arma::dot(tangent,x)*x;
        seg1_orth.col(lever_size-i) = tangent-seg1_para.col(lever_size-i);

        t_orthsq(lever_size-i) = arma::dot(seg1_orth.col(lever_size-i),seg1_orth.col(lever_size-i));
    }

    for (int i=0;i<lever_size;i++) {
        id_from = pmod(rem_hv2+i,num_bp);
        id_to   = pmod(rem_hvp+i,num_bp);
        triads->slice(id_to) = rot_align2*triads->slice(id_from);

        tangent = triads->slice(id_to).col(2)*disc_len;
        seg2_para.col(i) = arma::dot(tangent,x)*x;
        seg2_orth.col(i) = tangent-seg2_para.col(i);

        t_orthsq(i+lever_size) = arma::dot(seg2_orth.col(i),seg2_orth.col(i));
    }

//    for (int i=0;i<lever_size;i++) {
//        t_orthsq(i)            = arma::dot(seg1_orth.col(i),seg1_orth.col(i));
//        t_orthsq(i+lever_size) = arma::dot(seg2_orth.col(i),seg2_orth.col(i));
//    }

    /*
        Newton-Raphson
    */

    int twolever       = 2*lever_size;
    double disc_len_sq = disc_len*disc_len;
    double sqrt_fac;
    double num_sum;
    double denum_sum;

    arma::colvec v_parallel = arma::dot(x,rot_align1*rem_v1+rot_align2*rem_v2)*x;
    double v_totsq = twolever*disc_len;
    v_totsq = v_totsq*v_totsq;

    double a = std::sqrt( (v_totsq-len_rem_w*len_rem_w)/(v_totsq- arma::dot(v_parallel,v_parallel) ));
    double a_last = a;
    double da = 1;

//    cout << "#######################" << endl;
//    cout << "Starting Newton Raphson" << endl;

    int step = 0;
    while (std::abs(da)>REM_NEWTON_RAPHSON_EPS && step<REM_NEWTON_RAPHSON_MAX_STEPS) {
        step++;
        num_sum=-len_rem_w;
        denum_sum=0;
        for (int i=0;i<twolever;i++) {
            sqrt_fac = std::sqrt(disc_len_sq-a*a*t_orthsq(i));
            num_sum += sqrt_fac;
            denum_sum += a*t_orthsq(i)/sqrt_fac;
        }
        a += num_sum/denum_sum;
//        cout.precision(17);
//        cout << "#" << step << " a = " << a << endl;
//        cout.precision(17);
//        cout << "#" << step << " da = " << da << endl;
        da =a-a_last;
        a_last = a;
    }

    if (a!=a || a<0 || std::abs(da) > 1e-10) {
        return false;
    }

    /*
        rotate the individual triads
    */
    arma::colvec para,orth;
    arma::colvec n_from,n_to,rotax;
    arma::mat    rot_single;
    for (int i=0;i<lever_size;i++) {
        // first segment
        id_to   = pmod(rem_hAp+i,num_bp);

        orth = a*seg1_orth.col(i);
        para = std::sqrt(disc_len_sq-a*a*t_orthsq(i))*x;

        tangent = orth + para;
        n_from = triads->slice(id_to).col(2)/arma::norm(triads->slice(id_to).col(2));
        n_to   = tangent/arma::norm(tangent);
        rotax  = arma::cross(n_from,n_to);
        theta = std::asin(arma::norm(rotax));
        if (arma::dot(n_from,n_to)<0) {
            theta = M_PI-theta;
        }
        rot_single = getRotMat(theta*rotax/arma::norm(rotax));
        triads->slice(id_to) = rot_single*triads->slice(id_to);
    }
    for (int i=0;i<lever_size;i++) {
        // second segment
        id_to   = pmod(rem_hvp+i,num_bp);

        orth = a*seg2_orth.col(i);
        para = std::sqrt(disc_len_sq-a*a*t_orthsq(i+lever_size))*x;

        tangent = orth + para;
        n_from = triads->slice(id_to).col(2)/arma::norm(triads->slice(id_to).col(2));
        n_to   = tangent/arma::norm(tangent);
        rotax  = arma::cross(n_from,n_to);
        theta = std::asin(arma::norm(rotax));
        if (arma::dot(n_from,n_to)<0) {
            theta = M_PI-theta;
        }
        rot_single = getRotMat(theta*rotax/arma::norm(rotax));
        triads->slice(id_to) = rot_single*triads->slice(id_to);
    }
    triads->slice(rem_hBp) = triads->slice(rem_hB);

    /*
        reconstruct the positions
    */
    pos->col(rem_hAp) = pos->col(rem_hA);
    for (int i=1;i<twolever+1;i++) {
        id_to   = pmod(rem_hAp+i,num_bp);
        id_from = pmod(rem_hAp+i-1,num_bp);

        pos->col(id_to) = pos->col(id_from) + triads->slice(id_from).col(2)*disc_len;
    }
    // check if endpositions coincide
    id_from = pmod(rem_hBp-1,num_bp);
    if ( arma::norm((pos->col(id_from) + triads->slice(id_from).col(2)*disc_len)-pos->col(rem_hBp)) > REM_DISPLACEMENT_EPS) {
        cout << "REM endpoint discrepancy exceeds maximum (" << arma::norm((pos->col(id_from) + triads->slice(id_from).col(2)*disc_len)-pos->col(rem_hBp)) << ")" << endl;
        return false;
    }
    if (a<0.25) {
        return false;
    }
//    cout << "a = " << a << endl;
    return true;

}


void MCS_RepTherm::MID_move() {
/*
    translates the index of the segment inbetween the two modification regions.
*/
    int id1,id2;
    if (rep_forward) {
        for (int i=1;i<intseg_size;i++) {
            id1 = pmod(rem_hB+i,num_bp);
            id2 = pmod(rem_hBp+i,num_bp);

            pos->col(id2)       = pos->col(id1);
            triads->slice(id2)  = triads->slice(id1);
        }
    }
    else {
        for (int i=1;i<intseg_size;i++) {
            id1 = pmod(rem_hA-i,num_bp);
            id2 = pmod(rem_hAp-i,num_bp);

            pos->col(id2)       = pos->col(id1);
            triads->slice(id2)  = triads->slice(id1);
        }
    }
}

bool MCS_RepTherm::ADD_move() {

    /*
        Generate new segment
    */
    arma::mat  bp_pos_newseg(3,rep_size,arma::fill::zeros);
    arma::cube triads_newseg(3,3,rep_size);
    arma::colvec seg_vec;
    double       len_dw,len_dw_sq;


    len_dw  = rep_size*disc_len;
    len_dw_sq = len_dw*len_dw;

    /*
        Open ADD segment
    */

    arma::colvec w,v1,v2;
    arma::colvec nv1,nv2;
    arma::colvec x,y,z;

    double theta1,theta2,theta;
    double len_w,len_v1,len_v2;

    w  = pos->col(add_hB)-pos->col(add_hA);
    v1 = pos->col(add_hv)-pos->col(add_hA);
    v2 = pos->col(add_hB)-pos->col(add_hv);

    len_w  = arma::norm(w);
    len_v1 = arma::norm(v1);
    len_v2 = arma::norm(v2);

    nv1 = v1/len_v1;
    nv2 = v2/len_v2;

    x = w/len_w;
    y = v1-arma::dot(x,v1)*x;
    y = y/arma::norm(y);
    z = arma::cross(x,y);

    theta1 = std::acos(arma::dot(x,nv1));
    theta2 = std::acos(arma::dot(x,nv2));

    /*
        We find theta with Newton Raphson
    */

    /////////////////////////////////////////
    // Initial theta (theta_0)
    double theta_0, dtheta, theta_last;
    double alpha1,alpha2;
    double A,B,C,D;
    double cos_alpha1,cos_alpha2,sin_alpha1,sin_alpha2;
    double func, deriv;
    int step;

    /////////////////////////////////////////
    /////////////////////////////////////////
    theta_0 = add_NR_mean_theta0/add_NR_count_theta0;

    theta      = theta_0;
    theta_last = theta_0;
    dtheta = 4;

    step = 0;
    while (std::abs(dtheta)>ADD_NEWTON_RAPHSON_THETA_EPS && step<ADD_NEWTON_RAPHSON_MAX_STEPS) {

        /////////////////
        alpha1 = theta1+theta;
        alpha2 = theta2+theta;

        cos_alpha1 = std::cos(alpha1);
        sin_alpha1 = std::sin(alpha1);
        cos_alpha2 = std::cos(alpha2);
        sin_alpha2 = std::sin(alpha2);

        A = len_v1*cos_alpha1 + len_v2*cos_alpha2 - len_w;
        B = len_v1*sin_alpha1 - len_v2*sin_alpha2;
        C = len_v1*sin_alpha1 + len_v2*sin_alpha2;
        D = len_v2*cos_alpha2 - len_v1*cos_alpha1;

        func  = A*A+B*B-len_dw_sq;
        deriv = 2*(D*B-C*A);

        //////////////////////////
        // Newton Raphson Step ///
        theta      = theta - (func/deriv);
        //////////////////////////
        dtheta     = theta-theta_last;
        theta_last = theta;

        step++;
    }

    if (theta!=theta || std::abs(dtheta) > 1e-14) {
        return false;
    }

    /*
        Check Validity of solution and calculate rotation matrices
    */

    arma::colvec r1,r2;
    arma::colvec connect_dir;
    arma::mat    R1,R2;
    double       dist;

    R1 = getRotMat(theta*z);
    R2 = getRotMat(-theta*z);
    r1 = pos->col(add_hA) + R1*v1;
    r2 = pos->col(add_hB) - R2*v2;
    connect_dir = r2-r1;
    dist = arma::norm(connect_dir);
    connect_dir = connect_dir/dist;

    if (std::abs(dist-len_dw)>ADD_NEWTON_RAPHSON_DISPLACEMENT_EPS) {
        cout << "ADD point displacement distance exceeds maximum (" << std::abs(dist-len_dw) << ")" << endl;
        return false;
    }

    /*
        Calibrate theta_0 for NR
    */
    add_NR_mean_theta0 += theta;
    add_NR_count_theta0++;

    /*
        Rotate triads and positions from hinge add_hAp to add_hv1p
    */

    int id;     // id before
    int idp;    // id after (primed id)
    int idp_p1; // id after plus 1 (primed id+1)

    pos->col(add_hAp) = pos->col(add_hA);
    for (int i=0;i<lever_size;i++) {
        id      = pmod(add_hA  +i ,num_bp);
        idp     = pmod(add_hAp +i ,num_bp);
        idp_p1  = pmod(idp     +1 ,num_bp);

        triads->slice(idp) = R1*triads->slice(id);
        pos->col(idp_p1)   = pos->col(idp) + triads->slice(idp).col(2)*disc_len;
    }

    // Check consistency
    if (arma::norm(r1-pos->col(add_hv1p))>ADD_NEWTON_RAPHSON_DISPLACEMENT_EPS) {
        cout << "ADD point displacement distance of individual rotation of the first segment exceeds maximum (" << arma::norm(r1-pos->col(add_hv1p)) << ") " << endl;
        return false;
    }

    /*
        Rotate triads and positions from hinge add_hv2p to add_hBp
    */

    triads->slice(add_hBp)  = triads->slice(add_hB);
    pos->col(add_hBp)       = pos->col(add_hB);

    for (int i=1;i<=lever_size;i++) {
        id     = pmod(add_hB  -i ,num_bp);
        idp    = pmod(add_hBp -i ,num_bp);
        idp_p1 = pmod(idp     +1 ,num_bp);

        triads->slice(idp) = R2*triads->slice(id);
        pos->col(idp)   = pos->col(idp_p1) - triads->slice(idp).col(2)*disc_len;
    }

    if (arma::norm(r2-pos->col(add_hv2p))>ADD_NEWTON_RAPHSON_DISPLACEMENT_EPS) {
        cout << "ADD point displacement distance of individual rotation of the second segment exceeds maximum (" << arma::norm(r2-pos->col(add_hv2p)) << ") " << endl;
        return false;
    }

    /*
        Insert straight piece
    */

    // initial twist difference
//    double init_delta_twist = BPS[pmod(add_hv1p-1,num_bp)]->get_Theta_component(2);



    arma::colvec rotax;
    arma::mat    RotMat;
    double phi,len_rotax;
    int add_hv1p_m1 = pmod(add_hv1p-1,num_bp);

    rotax = arma::cross(triads->slice(add_hv1p_m1).col(2),connect_dir);
    len_rotax = arma::norm(rotax);
    rotax = rotax/len_rotax;
    phi   = std::asin(len_rotax);
    if (arma::dot(triads->slice(add_hv1p_m1).col(2),connect_dir)<0) {
        phi = M_PI-phi;
    }
    RotMat = getRotMat(phi*rotax);

    /*
        First rotate the triad at position add_hv1p -1 around the local tangent
        by the amount of hte intrinsic twist and then rotate it onto the vector
        connect_dir that connects the two open ends.
    */
    triads->slice(add_hv1p) =  RotMat * triads->slice(add_hv1p_m1) * Rotz(BPS[add_hv1p_m1]->get_T0_component(2));



    int id_prev;
    for (int i=add_hv1p+1;i<add_hv1p+rep_size;i++) {
        id_prev = pmod(i-1,num_bp);
        id      = pmod(i,num_bp);
        pos->col(id)      = pos->col(id_prev) + connect_dir*disc_len;
        triads->slice(id) = triads->slice(id_prev) * Rotz(BPS[id_prev]->get_T0_component(2));
    }
    return true;
}







bool MCS_RepTherm::ADD_move_FullNewtonRaphson() {
/*
    This one won't work without a full overhaul.
*/



    /*
        Generate new segment
    */

    arma::mat  bp_pos_newseg(3,rep_size,arma::fill::zeros);
    arma::cube triads_newseg(3,3,rep_size);
    arma::colvec seg_vec;
    double       len_dw;

//    arma::mat L  = chain->get_avg_chol();
//    arma::colvec Theta_gen;
//    int id;
//    for (int i=1;i<rep_size;i++) {
//        Theta_gen = {normaldist(gen),normaldist(gen),normaldist(gen)};
//        Theta_gen = L*Theta_gen;
//        id = pmod(add_hv1p+i-1,num_bp);
//
//        triads_newseg.slice(i) = triads_newseg.slice(i-1)* getRotMat(*BPS[id]->get_T0)  *getRotMat(Theta_gen);
//        bp_pos_newseg.col(i)   = bp_pos_newseg.col(i-1) + triads_newseg.slice(i-1)*disc_len;
//    }
//    seg_vec = bp_pos_newseg(rep_size-1);
//    seg_len = arma::norm(seg_vec);

    int id;
    for (int i=1;i<rep_size;i++) {
        id = pmod(add_hv1p+i-1,num_bp);

        triads_newseg.slice(i) = triads_newseg.slice(i-1)* getRotMat(*BPS[id]->get_T0());
        bp_pos_newseg.col(i)(2)   = i*disc_len;
    }
    seg_vec = bp_pos_newseg.col(rep_size-1);
    len_dw  = rep_size*disc_len;

    /*
        Open ADD segment
    */

    arma::colvec x;
    arma::colvec w;
    double       len_w;

    int twolever=2*lever_size;

    w     = pos->col(add_hB)-pos->col(add_hA);
    len_w = arma::norm(w);
    x     = w/len_w;


    arma::colvec tangent;
    arma::mat seg_para(3,twolever);
    arma::mat seg_orth(3,twolever);
    arma::vec seg_orth_sq(twolever);

    for (int i=0;i<twolever;i++) {
        id = pmod(add_hA+i,num_bp);

        seg_para.col(i) = arma::dot(triads->slice(id).col(2),x)*x*disc_len;
        seg_orth.col(i) = triads->slice(id).col(2)*disc_len-seg_para.col(i);
        seg_orth_sq(i)  = arma::dot(seg_orth.col(i),seg_orth.col(i));
    }

    /*
        Newton-Raphson
    */

    double disc_len_sq = disc_len*disc_len;
    double sqrt_fac;
    double num_sum;
    double denum_sum;

    double v_totsq = twolever*disc_len;
    v_totsq = v_totsq*v_totsq;
    double w_min_dw = len_w-len_dw;

    double a = std::sqrt( (v_totsq-(w_min_dw)*(w_min_dw))/(v_totsq- len_w*len_w ));
    double a_last = a;
    double da = 1;

    if (a!=a) {
        cout << (v_totsq-(w_min_dw)*(w_min_dw))/(v_totsq- len_w*len_w ) << endl;
        std::exit(0);
    }

    cout << "#######################" << endl;
    cout << "Starting Newton Raphson" << endl;
    cout << "#" << 0 << " a = " << a << endl;

    int step = 0;
    while (std::abs(da)>1e-15 && step<20) {
        step++;
        num_sum=-w_min_dw;
        denum_sum=0;
        for (int i=0;i<twolever;i++) {
            sqrt_fac = std::sqrt(disc_len_sq-a*a*seg_orth_sq(i));
            if (sqrt_fac!=sqrt_fac) {
                cout << a << endl;
                cout << disc_len_sq-a*a*seg_orth_sq(i) << endl;
            }

            num_sum += sqrt_fac;
            denum_sum += a*seg_orth_sq(i)/sqrt_fac;
        }
        a += num_sum/denum_sum;
        cout.precision(17);
        cout << "#" << step << " a = " << a << endl;
//        cout.precision(17);
//        cout << "#" << step << " da = " << da << endl;
        da =a-a_last;
        a_last = a;
    }

    return true;
}


bool MCS_RepTherm::Pre_Thermalize() {

    /*
        Thermalize the REM segment
    */

    int id_first,id_last,id_range;
    id_first = rem_hAp-THERMALIZE_DIST;
    id_range = 2*(lever_size+THERMALIZE_DIST);
    id_last  = id_first+id_range;

    std::uniform_int_distribution<> gen_idA{id_first, id_last-THERMALIZE_S_MIN};
    std::uniform_int_distribution<> gen_CS_size{THERMALIZE_S_MIN, lever_size};

    int therm_moves =THERMALIZE_MOVES_PER_RESIDUE*id_range;
    int idA,idB,CS_size;
    int moves = 0;

    while (moves<therm_moves) {
        CS_size = gen_CS_size(gen);
        idA     = gen_idA(gen);
        idB     = idA + CS_size;
        if (idB >id_last) {
            continue;
        }

//        if (idA<idThermFirst || idB > idThermLast) {
//            cout << "violation REM segment" << endl;
//            cout << idThermFirst << " " << idThermLast << endl;
//            cout << idA << " " << idB << endl;
//        }

//        cout << "Thermalize: " << pmod(idA,num_bp) << " " << pmod(idB,num_bp) << endl;
        MC_pre_thermal_rot->rotation(idA,idB);
        moves++;
    }

    /*
        Thermalize the ADD segment
    */
    id_first = add_hAp-THERMALIZE_DIST;
    id_range = 2*(lever_size+THERMALIZE_DIST)+rep_size;
    id_last  = id_first+id_range;

    decltype(gen_idA.param()) new_range_gen_idA(id_first, id_last-THERMALIZE_S_MIN);
    gen_idA.param(new_range_gen_idA);

    therm_moves = THERMALIZE_MOVES_PER_RESIDUE*id_range;
    moves = 0;

    while (moves<therm_moves) {
        CS_size = gen_CS_size(gen);
        idA     = gen_idA(gen);
        idB     = idA + CS_size;
        if (idB >id_last) {
            continue;
        }
//        cout << "Thermalize: " << pmod(idA,num_bp) << " " << pmod(idB,num_bp) << endl;
        MC_pre_thermal_rot->rotation(idA,idB);
        moves++;
    }

}


bool MCS_RepTherm::Thermalize() {

    /*
        Thermalize the REM segment
    */

    int id_first,id_last,id_range;
    id_first = rem_hAp-THERMALIZE_DIST;
    id_range = 2*(lever_size+THERMALIZE_DIST);
    id_last  = id_first+id_range;

    std::uniform_int_distribution<> gen_idA{id_first, id_last-THERMALIZE_S_MIN};
    std::uniform_int_distribution<> gen_CS_size{THERMALIZE_S_MIN, rep_size};

    int therm_moves =THERMALIZE_MOVES_PER_RESIDUE*id_range;
    int idA,idB,CS_size;
    int moves = 0;

    int accept_count = 0;

    while (moves<therm_moves) {
        CS_size = gen_CS_size(gen);
        idA     = gen_idA(gen);
        idB     = idA + CS_size;
        if (idB >id_last) {
            continue;
        }
        if (MC_thermal_rot->rotation(idA,idB)) {
            accept_count++;
        }
        moves++;
    }

    /*
        Thermalize the ADD segment
    */
    id_first = add_hAp-THERMALIZE_DIST;
    id_range = 2*(lever_size+THERMALIZE_DIST)+rep_size;
    id_last  = id_first+id_range;

    decltype(gen_idA.param()) new_range_gen_idA(id_first, id_last-THERMALIZE_S_MIN);
    gen_idA.param(new_range_gen_idA);

    therm_moves = THERMALIZE_MOVES_PER_RESIDUE*id_range;
    moves = 0;

    while (moves<therm_moves) {
        CS_size = gen_CS_size(gen);
        idA     = gen_idA(gen);
        idB     = idA + CS_size;
        if (idB >id_last) {
            continue;
        }
        if (MC_thermal_rot->rotation(idA,idB)) {
            accept_count++;
        }
        moves++;
    }
//    cout << "Acceptance Rate: " << 1.*accept_count/moves << endl;

}


bool MCS_RepTherm::Thermalize_Hinges() {

    /*
        Thermalize remAp hinge
    */

    int id_first,id_last,id_range;
    id_first = rem_hAp-THERMALIZE_DIST;
    id_range = 2*(THERMALIZE_DIST);
    id_last  = id_first+id_range;

    std::uniform_int_distribution<> gen_idA{id_first, id_last-THERMALIZE_S_MIN};
    std::uniform_int_distribution<> gen_CS_size{THERMALIZE_S_MIN, THERMALIZE_DIST};

    int therm_moves =THERMALIZE_MOVES_PER_RESIDUE*id_range;
    int idA,idB,CS_size;
    int moves = 0;

    int accept_count = 0;

    while (moves<therm_moves) {
        CS_size = gen_CS_size(gen);
        idA     = gen_idA(gen);
        idB     = idA + CS_size;
        if (idB >id_last) {
            continue;
        }
        if (MC_thermal_rot->rotation(idA,idB)) {
            accept_count++;
        }
        moves++;
    }

    /*
        Thermalize remvp hinge
    */
    id_first = rem_hvp-THERMALIZE_DIST;
    id_range = 2*(THERMALIZE_DIST);
    id_last  = id_first+id_range;

    decltype(gen_idA.param()) new_range_gen_remvp(id_first, id_last-THERMALIZE_S_MIN);
    gen_idA.param(new_range_gen_remvp);

    therm_moves = THERMALIZE_MOVES_PER_RESIDUE*id_range;


    while (moves<therm_moves) {
        CS_size = gen_CS_size(gen);
        idA     = gen_idA(gen);
        idB     = idA + CS_size;
        if (idB >id_last) {
            continue;
        }
        if (MC_thermal_rot->rotation(idA,idB)) {
            accept_count++;
        }
        moves++;
    }

    /*
        Thermalize remBp hinge
    */
    id_first = rem_hBp-THERMALIZE_DIST;
    id_range = 2*(THERMALIZE_DIST);
    id_last  = id_first+id_range;

    decltype(gen_idA.param()) new_range_gen_remBp(id_first, id_last-THERMALIZE_S_MIN);
    gen_idA.param(new_range_gen_remBp);

    therm_moves = THERMALIZE_MOVES_PER_RESIDUE*id_range;


    while (moves<therm_moves) {
        CS_size = gen_CS_size(gen);
        idA     = gen_idA(gen);
        idB     = idA + CS_size;
        if (idB >id_last) {
            continue;
        }
        if (MC_thermal_rot->rotation(idA,idB)) {
            accept_count++;
        }
        moves++;
    }

    /*
        Thermalize addAp hinge
    */
    id_first = add_hAp-THERMALIZE_DIST;
    id_range = 2*(THERMALIZE_DIST);
    id_last  = id_first+id_range;

    decltype(gen_idA.param()) new_range_gen_addAp(id_first, id_last-THERMALIZE_S_MIN);
    gen_idA.param(new_range_gen_addAp);

    therm_moves = THERMALIZE_MOVES_PER_RESIDUE*id_range;


    while (moves<therm_moves) {
        CS_size = gen_CS_size(gen);
        idA     = gen_idA(gen);
        idB     = idA + CS_size;
        if (idB >id_last) {
            continue;
        }
        if (MC_thermal_rot->rotation(idA,idB)) {
            accept_count++;
        }
        moves++;
    }

    /*
        Thermalize added segment around hinges addv1p and addv2p
    */
    id_first = add_hv1p-THERMALIZE_DIST;
    id_range = 2*(THERMALIZE_DIST)+rep_size;
    id_last  = id_first+id_range;

    decltype(gen_idA.param()) new_range_gen_addedseg(id_first, id_last-THERMALIZE_S_MIN);
    gen_idA.param(new_range_gen_addedseg);

    therm_moves = THERMALIZE_MOVES_PER_RESIDUE*id_range;


    while (moves<therm_moves) {
        CS_size = gen_CS_size(gen);
        idA     = gen_idA(gen);
        idB     = idA + CS_size;
        if (idB >id_last) {
            continue;
        }
        if (MC_thermal_rot->rotation(idA,idB)) {
            accept_count++;
        }
        moves++;
    }

    /*
        Thermalize addBp hinge
    */
    id_first = add_hBp-THERMALIZE_DIST;
    id_range = 2*(THERMALIZE_DIST);
    id_last  = id_first+id_range;

    decltype(gen_idA.param()) new_range_gen_addBp(id_first, id_last-THERMALIZE_S_MIN);
    gen_idA.param(new_range_gen_addBp);

    therm_moves = THERMALIZE_MOVES_PER_RESIDUE*id_range;

    while (moves<therm_moves) {
        CS_size = gen_CS_size(gen);
        idA     = gen_idA(gen);
        idB     = idA + CS_size;
        if (idB >id_last) {
            continue;
        }
        if (MC_thermal_rot->rotation(idA,idB)) {
            accept_count++;
        }
        moves++;
    }

//    cout << "Acceptance Rate: " << 1.*accept_count/moves << endl;

}



bool MCS_RepTherm::Check_EV() {

    /*
        REM segment
    */
    int rem_first,rem_last;
    int add_first,add_last;

    rem_first = pmod(rem_hAp-THERMALIZE_DIST,num_bp);
    rem_last  = pmod(rem_first+2*(lever_size+THERMALIZE_DIST),num_bp);

    /*
        ADD segment
    */
    add_first = pmod(add_hAp-THERMALIZE_DIST,num_bp);
    add_last  = pmod(add_first+2*(lever_size+THERMALIZE_DIST)+rep_size,num_bp);


    int B1,B2,D1,D2;
    if (rep_forward) {
        B1 = rem_first;
        B2 = rem_last;
        D1 = add_first;
        D2 = add_last;
    }
    else {
        B1 = add_first;
        B2 = add_last;
        D1 = rem_first;
        D2 = rem_last;
    }

    /*
        Situations I,II and III

    */
    if (B1 < D2) {

        moved_intervals[0](EV_FROM) = 0;
        moved_intervals[0](EV_TO)   = B1-1;
        if (B1 >0)  moved_intervals[0](EV_TYPE) = 0;
        else        moved_intervals[0](EV_TYPE) = -1;

        moved_intervals[1](EV_FROM) = B1;
        moved_intervals[1](EV_TO)   = B2;
        moved_intervals[1](EV_TYPE) = 2100;

        moved_intervals[2](EV_FROM) = B2+1;
        moved_intervals[2](EV_TO)   = D1-1;
        moved_intervals[2](EV_TYPE) = 2;

        moved_intervals[3](EV_FROM) = D1;
        moved_intervals[3](EV_TO)   = D2;
        moved_intervals[3](EV_TYPE) = 2101;

        moved_intervals[4](EV_FROM) = D2+1;
        moved_intervals[4](EV_TO)   = num_bp-1;
        if (D2 == num_bp-1) moved_intervals[4](EV_TYPE) = -1;
        else                moved_intervals[4](EV_TYPE) = 0;

    }

    else {
        /*
            Cases where neither B nore D are split (Situation IV,V and VI)
        */
        if (B1<B2 && D1<D2){
            moved_intervals[0](EV_FROM) = 0;
            moved_intervals[0](EV_TO)   = D1-1;
            if (D1 > 0) moved_intervals[0](EV_TYPE) = 2;
            else        moved_intervals[0](EV_TYPE) = -1;

            moved_intervals[1](EV_FROM) = D1;
            moved_intervals[1](EV_TO)   = D2;
            moved_intervals[1](EV_TYPE) = 2100;

            moved_intervals[2](EV_FROM) = D2+1;
            moved_intervals[2](EV_TO)   = B1-1;
            moved_intervals[2](EV_TYPE) = 0;

            moved_intervals[3](EV_FROM) = B1;
            moved_intervals[3](EV_TO)   = B2;
            moved_intervals[3](EV_TYPE) = 2101;

            moved_intervals[4](EV_FROM) = B2+1;
            moved_intervals[4](EV_TO)   = num_bp-1;
            if (B2 == num_bp-1) moved_intervals[4](EV_TYPE) = -1;
            else                moved_intervals[4](EV_TYPE) = 2;
        }
        else {
            /*
                D is split (sitatuation VII)
            */
            if (D1>D2) {
                moved_intervals[0](EV_FROM) = 0;
                moved_intervals[0](EV_TO)   = D2;
                moved_intervals[0](EV_TYPE) = 2100;

                moved_intervals[1](EV_FROM) = D2+1;
                moved_intervals[1](EV_TO)   = B1-1;
                moved_intervals[1](EV_TYPE) = 0;

                moved_intervals[2](EV_FROM) = B1;
                moved_intervals[2](EV_TO)   = B2;
                moved_intervals[2](EV_TYPE) = 2101;

                moved_intervals[3](EV_FROM) = B2+1;
                moved_intervals[3](EV_TO)   = D1-1;
                moved_intervals[3](EV_TYPE) = 2;

                moved_intervals[4](EV_FROM) = D1;
                moved_intervals[4](EV_TO)   = num_bp-1;
                moved_intervals[4](EV_TYPE) = 2100;
            }
            else {
                /*
                    B is split (sitatuation VIII)
                */
                if (B1>B2) {
                    moved_intervals[0](EV_FROM) = 0;
                    moved_intervals[0](EV_TO)   = B2;
                    moved_intervals[0](EV_TYPE) = 2100;

                    moved_intervals[1](EV_FROM) = B2+1;
                    moved_intervals[1](EV_TO)   = D1-1;
                    moved_intervals[1](EV_TYPE) = 2;

                    moved_intervals[2](EV_FROM) = D1;
                    moved_intervals[2](EV_TO)   = D2;
                    moved_intervals[2](EV_TYPE) = 2101;

                    moved_intervals[3](EV_FROM) = D2+1;
                    moved_intervals[3](EV_TO)   = B1-1;
                    moved_intervals[3](EV_TYPE) = 0;

                    moved_intervals[4](EV_FROM) = B1;
                    moved_intervals[4](EV_TO)   = num_bp-1;
                    moved_intervals[4](EV_TYPE) = 2100;
                }
                else {
                    cout << "I missed something!" << endl;
                    std::exit(0);
                }
            }
        }

    }


    EV->set_backup_conf(&backup_bp_pos1,&backup_triads1);

    return EV->check(&moved_intervals);
}






