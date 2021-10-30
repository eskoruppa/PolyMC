#include "Dump_HatCurve.h"

Dump_HatCurve::Dump_HatCurve(Chain * ch, int N_dump, const std::string& filename, int N_dump_to_file, bool append)
: Dump(ch,N_dump,filename,append),
N_dump_to_file(N_dump_to_file/N_dump)
{
    z       = 0;
    zsq     = 0;
    x       = 0;
    xsq     = 0;
    y       = 0;
    ysq     = 0;
    Wr      = 0;
    Wrsq    = 0;
    counter = 0;

}
Dump_HatCurve::~Dump_HatCurve() {}


void Dump_HatCurve::prod_dump() {
    double new_z  = arma::dot(chain->get_force_dir(),chain->get_bp_pos()->col(num_bp-1)-chain->get_bp_pos()->col(0));
    double new_x  = chain->get_bp_pos()->col(num_bp-1)(0)-chain->get_bp_pos()->col(0)(0);
    double new_y  = chain->get_bp_pos()->col(num_bp-1)(1)-chain->get_bp_pos()->col(0)(1);
    double new_Tw = chain->cal_twist(0,num_bp);
    double new_Wr = chain->get_dLK() - new_Tw;
    z    += new_z;
    zsq  += new_z*new_z;
    x    += new_x;
    xsq  += new_x*new_x;
    y    += new_y;
    ysq  += new_y*new_y;
    Wr   += new_Wr;
    Wrsq += new_Wr*new_Wr;
    Tw   += new_Tw;
    Twsq += new_Tw*new_Tw;

    counter++;

    if (counter >=N_dump_to_file) {
        std::ofstream ofstr;
        ofstr.open(fn, std::ofstream::out | std::ofstream::app);

//        double cal_Wr = chain->cal_langowski_writhe_1a(1);
//        double cal_Tw = chain->cal_twist(0,num_bp);
//        double cal_Lk = cal_Wr + cal_Tw;
        double dLK_fix;
        if (chain->torsional_trap_active()) {
            dLK_fix = chain->get_torstrap_dLK_fix();
        }
        else {
            dLK_fix = chain->get_dLK();
        }

        ofstr << chain->get_force() << " " << dLK_fix << " " << chain->get_disc_len()*num_bps << " " << z/counter << " " << zsq/counter << " " << x/counter << " " << xsq/counter << " " << y/counter << " " << ysq/counter << " " << Wr/counter  << " " << Wrsq/counter << " " << Tw/counter  << " " << Twsq/counter  << std::endl;
        ofstr.close();
        z       = 0;
        zsq     = 0;
        x       = 0;
        xsq     = 0;
        y       = 0;
        ysq     = 0;
        Wr      = 0;
        Wrsq    = 0;
        Tw      = 0;
        Twsq    = 0;
        counter = 0;
    }
}

void Dump_HatCurve::final_dump() {
}



////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////




Dump_HatCurveStatistics::Dump_HatCurveStatistics(Chain * ch, int N_dump, std::string filename, bool append)
: Dump(ch,N_dump,filename,append) {}
Dump_HatCurveStatistics::~Dump_HatCurveStatistics() {}


void Dump_HatCurveStatistics::prod_dump() {
    double new_z  = arma::dot(chain->get_force_dir(),chain->get_bp_pos()->col(num_bp-1)-chain->get_bp_pos()->col(0));
    double new_Wr = chain->get_dLK() - chain->cal_twist(0,num_bp);

    std::ofstream ofstr;
    ofstr.open(fn, std::ofstream::out | std::ofstream::app);

    ofstr <<  new_z << " " << new_Wr << std::endl;
    ofstr.close();
}

void Dump_HatCurveStatistics::final_dump() {
}




