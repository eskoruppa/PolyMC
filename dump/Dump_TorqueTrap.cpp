#include "Dump_TorqueTrap.h"

Dump_TorqueTrap::Dump_TorqueTrap(Chain * ch, int N_dump, const std::string& filename, bool append)
: Dump(ch,N_dump,filename,append)
{
    if (!app) {
        std::ofstream ofstr;
        ofstr.open(fn, std::ofstream::out | std::ofstream::trunc);

        double trapstiff, trap_dLK_fix, L, f;

        trapstiff    = chain->get_torstrap_trapstiff();
        trap_dLK_fix = chain->get_torstrap_dLK_fix();
        L            = chain->get_disc_len()*num_bps;
        f            = chain->get_force();

        ofstr << f << " " << L << " " << trap_dLK_fix << " " << trapstiff << std::endl;
        ofstr.close();
    }
}
Dump_TorqueTrap::~Dump_TorqueTrap() {}


void Dump_TorqueTrap::prod_dump() {

    double torque, Lk, Tw, zext;
//    torque = chain->measure_current_torque();
    torque = chain->get_torque_measure_accu(true);
    Lk = chain->get_dLK();
    Tw = chain->cal_twist(0,num_bp);
    zext = arma::dot(chain->get_force_dir(),chain->get_bp_pos()->col(num_bp-1)-chain->get_bp_pos()->col(0));

	std::ofstream ofstr;
	ofstr.open(fn, std::ofstream::out | std::ofstream::app);
//    ofstr << chain->get_torque_measure_accu(true) << std::endl;
    ofstr << zext << " " << torque << " " << Lk << " " << Tw << std::endl;
    ofstr.close();

//	ofstream ofstr1;
//	ofstr1.open(fn, ofstream::out | ofstream::app);
//    ofstr1 << chain->cal_twist(0,num_bp) << endl;
//    ofstr1.close();
}

void Dump_TorqueTrap::final_dump() {
}

