#include "Dump_Endbead.h"

/////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////

Dump_Endbead::Dump_Endbead(Chain * ch, int N_dump, const std::string& filename, bool append)
: Dump(ch,N_dump,filename,append)
{
    if (!app) {
        std::ofstream ofstr;
        ofstr.open(fn, std::ofstream::out | std::ofstream::trunc);
        ofstr.close();
    }
    init_store(DUMP_ENDBEAD_DEFAULT_STORE_NUM);
}
Dump_Endbead::Dump_Endbead(Chain * ch, int N_dump, int N_to_file, const std::string& filename, bool append)
: Dump(ch,N_dump,filename,append)
{
    if (!app) {
        std::ofstream ofstr;
        ofstr.open(fn, std::ofstream::out | std::ofstream::trunc);
        ofstr.close();
    }
    init_store(N_to_file);
}
Dump_Endbead::~Dump_Endbead() {}

void Dump_Endbead::init_store(int num) {
    store_num     = num;
    store_topol   = arma::zeros(2,num);
    store_counter = 0;
}

void Dump_Endbead::write2file() {
	std::ofstream ofstr;
	ofstr.open(fn, std::ofstream::out | std::ofstream::app);
	for (unsigned i=0;i<store_counter;i++) {
        ofstr << store_topol(0,i) << " " << store_topol(1,i) << std::endl;
    }
    ofstr.close();
    store_counter = 0;
}

void Dump_Endbead::prod_dump() {
    store_topol(0,store_counter) = chain->get_dLK();
    store_topol(1,store_counter) = chain->cal_twist();
    store_counter++;
    if (store_counter == store_num) {
        write2file();
    }
}

void Dump_Endbead::final_dump() {
    write2file();
}


/////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////


Dump_Endbead_Stats::Dump_Endbead_Stats(Chain * ch, int N_dump, const std::string& filename, bool append)
: Dump(ch,1,filename,append)
{
    if (!app) {
        std::ofstream ofstr;
        ofstr.open(fn, std::ofstream::out | std::ofstream::trunc);
        ofstr << chain->get_force() << " " << chain->get_contour_len() << std::endl;
        ofstr.close();
    }

    sum_dump_every = N_dump;
    sum       = 0;
    sqsum     = 0;
    count_sum = 0;
}
Dump_Endbead_Stats::~Dump_Endbead_Stats() {}

void Dump_Endbead_Stats::prod_dump() {

    sum   += chain->get_dLK();
    sqsum += chain->get_dLK()*chain->get_dLK();
    count_sum ++;

    if (count_sum >= sum_dump_every) {
        std::ofstream ofstr;
        ofstr.open(fn, std::ofstream::out | std::ofstream::app);
        ofstr << sum/count_sum << " " << sqsum/count_sum << std::endl;
        ofstr.close();
        sum       = 0;
        sqsum     = 0;
        count_sum = 0;
    }
}

void Dump_Endbead_Stats::final_dump() {
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////


Dump_Endbead_Ceff::Dump_Endbead_Ceff(Chain * ch, int N_dump, const std::string& filename, bool append)
: Dump(ch,N_dump,filename,append)
{
    if (!app) {
        std::ofstream ofstr;
        ofstr.open(fn, std::ofstream::out | std::ofstream::trunc);
        ofstr.close();
    }
    sum       = 0;
    sqsum     = 0;
    count_sum = 0;
}
Dump_Endbead_Ceff::~Dump_Endbead_Ceff() {}

void Dump_Endbead_Ceff::prod_dump() {
    sum   += chain->get_dLK();
    sqsum += chain->get_dLK()*chain->get_dLK();
    count_sum ++;
}

void Dump_Endbead_Ceff::final_dump() {
    double var = sqsum/count_sum - sum*sum/count_sum/count_sum;
    double Ceff = chain->get_contour_len()/(4*M_PI*M_PI*var);

    std::ofstream ofstr;
    ofstr.open(fn, std::ofstream::out | std::ofstream::app);
    ofstr << Ceff << " " << chain->get_force() << " " << chain->get_contour_len() << " " << sum/count_sum << " " << sqsum/count_sum << std::endl;
    ofstr.close();

}


