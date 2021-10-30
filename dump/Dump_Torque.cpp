#include "Dump_Torque.h"

Dump_Torque::Dump_Torque(Chain * ch, int N_dump, const std::string& filename, bool append)
: Dump(ch,N_dump,filename,append)
{
    if (!app) {
        std::ofstream ofstr;
        ofstr.open(fn, std::ofstream::out | std::ofstream::trunc);
        ofstr.close();
    }
}
Dump_Torque::~Dump_Torque() {}


void Dump_Torque::prod_dump() {
	std::ofstream ofstr;
	ofstr.open(fn, std::ofstream::out | std::ofstream::app);
    ofstr << chain->get_torque_measure_accu(true) << std::endl;
    ofstr.close();

//	ofstream ofstr1;
//	ofstr1.open(fn, ofstream::out | ofstream::app);
//    ofstr1 << chain->cal_twist(0,num_bp) << endl;
//    ofstr1.close();
}

void Dump_Torque::final_dump() {
}

