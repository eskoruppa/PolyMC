#include "../Chain.h"


void Chain::allow_free_endbead_rotation() {
    fix_link(false);
    link_torsional_trap = false;
    link_track_terminus = true;
    dLK = cal_langowski_writhe_1a()+cal_twist(0,num_bps-1);
    terminus_ref_triad = triads.slice(num_bp-1);
    change_torque(0);
}

////////////////////// CURRENT LINKING NUMBER //////////////////////

double Chain::get_dLK() {
    return dLK;
}

arma::mat* Chain::get_terminus_ref_triad() {
    return &terminus_ref_triad;
}

void Chain::set_dLK(double dLK) {
    this->dLK = dLK;
}

void Chain::set_terminus_ref_triad(const arma::mat & ref_triad) {
    this->terminus_ref_triad = ref_triad;
}



////////////////////////////// OPTIONS //////////////////////////////

bool Chain::torsional_trap_active() {
    return link_torsional_trap;
}
bool Chain::const_torque_active() {
    return link_const_torque;
}
bool Chain::track_terminus_active() {
    return link_track_terminus;
}

void Chain::fix_link(bool fix) {
    link_constrained = fix;
}
bool Chain::link_fixed() {
    return link_constrained;
}

////////////////////////////// TORSIONAL TRAP //////////////////////////////

void Chain::impose_torsional_trap(double dLK_fix, double trapstiff) {

    dLK = cal_langowski_writhe_1a()+cal_twist(0,num_bps-1);
    terminus_ref_triad = triads.slice(num_bp-1);

    fix_link(false);
    link_track_terminus = true;
    link_const_torque   = false;
    if (trapstiff==0) { // || dLK_fix==0) {
        link_torsional_trap = false;
        std::cout << "Torsional Trap deactivated!" << std::endl;
    }
    else {
        std::cout << "Activated Torsional Trap" << std::endl;

        link_torsional_trap = true;
        torstrap_dLK_fix = dLK_fix;

        if (trapstiff>0) {
            torstrap_stiff = trapstiff;
            torstrap_dLK_aim = 15./360;
        }
        else {

            double meanC = 0;
            //arma::mat stiffmat;
            for (unsigned i=0;i<num_bps;i++) {
                //stiffmat = arma::inv(BPS[i]->get_cov())*disc_len;
                meanC += disc_len/BPS[i]->get_cov()(2,2);
            }
            meanC /= num_bps;
            std::cout << " Mean Renormalized Torsional Stiffness (from diagonal components) = " << meanC << std::endl;

            /*
                TODO: dLK_fix sets link_torsional_trap to false
            */
            if (torstrap_dLK_fix == 0) {
                torstrap_dLK_aim = 15./360;
                torstrap_stiff   = 4*M_PI*M_PI *meanC * kT_ref / get_contour_len()*100;
            }
            else {

                torstrap_dLK_aim = -sgn(torstrap_dLK_fix)*15./360.;
                if (std::abs(torstrap_dLK_aim/torstrap_dLK_fix) > 0.005) {
                    torstrap_dLK_aim = -torstrap_dLK_fix*0.005;
                }
                torstrap_stiff = 4*M_PI*M_PI *meanC * kT_ref / get_contour_len() * (std::abs(torstrap_dLK_fix)-std::abs(torstrap_dLK_aim))/(std::abs(torstrap_dLK_aim));
            }
            std::cout << " Linking Number fix = " << torstrap_dLK_fix << std::endl;
            std::cout << " Fluctuation aim    = " << torstrap_dLK_aim*360.0 << " deg" << std::endl;
        }
        torstrap_torque_sum   = 0;
        torstrap_torque_count = 0;
    }
}

double Chain::get_torstrap_dLK_fix() {
    return torstrap_dLK_fix;
}

void Chain::set_torstrap_dLK_fix(double dLK) {
    torstrap_dLK_fix = dLK;
}

double Chain::get_torstrap_dLK_aim() {
    return torstrap_dLK_aim;
}

double Chain::get_torstrap_trapstiff() {
    return torstrap_stiff;
}

double Chain::measure_current_torque() {
    return -torstrap_stiff*(dLK-torstrap_dLK_fix);
}

double Chain::get_torque_measure_accu(bool reset_sum) {
    if (torstrap_torque_count==0) {
        return 0;
    }
    double cal_tor = torstrap_torque_sum/torstrap_torque_count;
    if (reset_sum) {
        torstrap_torque_sum   = 0;
        torstrap_torque_count = 0;
    }
    return cal_tor;
}


////////////////////////////// TORQUE //////////////////////////////


void Chain::impose_torque(double torque) {

    fix_link(false);
    link_torsional_trap = false;
    link_track_terminus = true;

    dLK = cal_langowski_writhe_1a()+cal_twist(0,num_bps-1);
    terminus_ref_triad = triads.slice(num_bp-1);

    change_torque(torque);
}

void Chain::change_torque(double torque) {
    /*
        This method changes the torque independent on whether torque constraints are active or not.
    */
    if (torque==0) {
        link_const_torque = false;
        this->torque      = 0;
        this->beta_torque = 0;
    }
    else {
        link_const_torque = true;
        this->torque      = torque;
        this->beta_torque = torque/kT;
    }

}

double Chain::extract_torque_betaenergy() {
    if (link_const_torque) {
        return -2*M_PI*beta_torque * (dLK);
    }
    return 0;
}

/////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////  TERMINUS ROTATION  //////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////

double Chain::propose_terminus_rotation(const arma::mat & T_newlast, double change_angle) {
    /*
        Perhaps build in another boolian to indicate whether the terminus rotations has to be tracked
    */
    if (!link_track_terminus) {
        return 0;
    }

    #ifdef CHAIN_DEBUG
    if (std::abs(1-arma::dot(terminus_ref_triad.col(2),T_newlast.col(2))>1e-10)) {
        std::cout << terminus_ref_triad;
        std::cout << T_newlast;
        std::cout << "Error: Tracing terminus rotation for non-preserved termini!" << std::endl;
        std::exit(0);
    }
    #endif


//    dLK_propose = ExtractTheta3(terminus_ref_triad.t() * T_newlast )*ONE_OVER_TWO_PI+torstrap_dLK_fix;
//    dLK_propose = dLK_propose - std::round(dLK_propose - dLK);

    double dLK_check = dLK + change_angle*ONE_OVER_TWO_PI;
    dLK_propose = ExtractTheta3(terminus_ref_triad.t() * T_newlast )*ONE_OVER_TWO_PI;
    dLK_propose = dLK_propose + std::round(dLK_check - dLK_propose);

    /*
    ---------------------------------------
    */
    if (dLK_propose!=dLK_propose) {
        return 1e10;
        std::exit(0);
    }

    if (std::abs(dLK_propose-dLK) > M_PI) {
        return 1e10;
    }
    /*
    ---------------------------------------
    */


    /*
        Torsional Trap
    */
    if (link_torsional_trap) {
        double old_del_dLK,new_del_dLK;
        old_del_dLK = dLK         - torstrap_dLK_fix;
        new_del_dLK = dLK_propose - torstrap_dLK_fix;
        return beta*torstrap_stiff*0.5*(new_del_dLK*new_del_dLK - old_del_dLK*old_del_dLK);
    }

    /*
        Torque Energy
    */
    // TODO define constant for 2pi beta torque
    if (link_const_torque) {
        return -2*M_PI*beta_torque * (dLK_propose-dLK);
    }
    return 0;
}

void   Chain::set_terminus_rotation(bool accept) {
    if (accept) {
        dLK = dLK_propose;
    }
    if (link_torsional_trap) {
        torstrap_torque_sum  += measure_current_torque();
        torstrap_torque_count++;
    }
}

