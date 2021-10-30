#include "Dump_ForceExtension.h"

Dump_ForceExtension::Dump_ForceExtension(Chain * ch, int N_dump, const std::string& filename, int start_id, int end_id, bool append)
: Dump(ch,N_dump,filename,append)
{
    z      = 0;
    z_sq   = 0;
    counter = 0;
    iID = start_id;
    if (iID < 0) iID = 0;
    fID = end_id;
    if (fID > num_bp-1 || fID <= iID) fID = num_bp-1;
}
Dump_ForceExtension::~Dump_ForceExtension() {}


void Dump_ForceExtension::prod_dump() {
    double new_z = arma::dot(chain->get_force_dir(),chain->get_bp_pos()->col(fID)-chain->get_bp_pos()->col(iID));
    z    += new_z;
    z_sq += new_z*new_z;
    counter++;
}


void Dump_ForceExtension::final_dump() {
    z    /= counter;
    z_sq /= counter;
    double L = chain->get_disc_len()*(fID-iID);

	std::ofstream ofstr;
	if (app) ofstr.open(fn, std::ofstream::out | std::ofstream::app);
	else     ofstr.open(fn, std::ofstream::out | std::ofstream::trunc);

    ofstr << chain->get_force() << " " << counter << " " << z << " " << z_sq << " " << L << std::endl;
    ofstr.close();
}

