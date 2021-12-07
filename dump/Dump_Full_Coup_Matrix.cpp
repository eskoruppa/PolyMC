#include "Dump_Full_Coup_Matrix.h"

Dump_Full_Coup_Matrix::Dump_Full_Coup_Matrix(Chain * ch, const std::string& filename)
: Dump(ch,0,filename,false)
{
    int dim = num_bps*3;
    arma::mat full_mat = arma::zeros(dim,dim);
    std::vector<double> params;
    arma::mat part_mat;
    for (int i=0;i<num_bps;i++) {
        params = ch->get_BPS()->at(i)->get_eval_energy_params(0);
        full_mat(arma::span(i*3,(i+1)*3 - 1),arma::span(i*3,(i+1)*3 - 1)) = params2mat(params);

        for (int nl=1;i-nl>0;nl++) {
            params = ch->get_BPS()->at(i)->get_eval_energy_params(-nl);
            if (params.size() <= 0) {
                continue;
            }
            part_mat = params2mat(params);

            if (arma::accu(full_mat(arma::span((i-nl)*3,(i-nl+1)*3 - 1),arma::span(i*3,(i+1)*3 - 1)) - part_mat) > 1e-10) {
                std::cout << "Matrices inconsistent" << std::endl;
                std::exit(0);
            }
            if (arma::accu(full_mat(arma::span(i*3,(i+1)*3 - 1),arma::span((i-nl)*3,(i-nl+1)*3 - 1)) - part_mat.t()) > 1e-10) {
                std::cout << "Matrices inconsistent" << std::endl;
                std::exit(0);
            }
//            full_mat(arma::span((i-nl)*3,(i-nl+1)*3 - 1),arma::span(i*3,(i+1)*3 - 1)) = part_mat;
//            full_mat(arma::span(i*3,(i+1)*3 - 1),arma::span((i-nl)*3,(i-nl+1)*3 - 1)) = part_mat.t();
        }

        for (int nr=1;i+nr<num_bps;nr++) {
            params = ch->get_BPS()->at(i)->get_eval_energy_params(nr);
            if (params.size() <= 0) {
                continue;
            }
            part_mat = params2mat(params);
            full_mat(arma::span(i*3,(i+1)*3 - 1),arma::span((i+nr)*3,(i+nr+1)*3 - 1)) = part_mat;
            full_mat(arma::span((i+nr)*3,(i+nr+1)*3 - 1),arma::span(i*3,(i+1)*3 - 1)) = part_mat.t();
        }
    }
	std::ofstream ofstr;
	ofstr.open(fn+".idb_mat", std::ofstream::out | std::ofstream::trunc);

//	int row,col;

	for (int r=0;r<dim;r++) {
        for (int c=0;c<dim-1;c++) {
            ofstr << full_mat(r,c) <<  " ";
        }
        ofstr << full_mat(r,dim-1) << std::endl;
	}
    ofstr.close();


	ofstr.open(fn+".idb_T0", std::ofstream::out | std::ofstream::trunc);
    arma::colvec T0;
    for (int i=0;i<num_bps;i++) {
        T0 = *ch->get_BPS()->at(i)->get_T0();
        ofstr << i << " " << T0(0) <<  " " << T0(1) <<  " " << T0(2) << std::endl;
    }
    ofstr.close();




}


Dump_Full_Coup_Matrix::~Dump_Full_Coup_Matrix() {}


arma::mat Dump_Full_Coup_Matrix::params2mat(std::vector<double> params) {
    arma::mat mat = arma::zeros(3,3);
    mat(0,0) = params[0];
    mat(0,1) = params[1];
    mat(0,2) = params[2];
    mat(1,0) = params[3];
    mat(1,1) = params[4];
    mat(1,2) = params[5];
    mat(2,0) = params[6];
    mat(2,1) = params[7];
    mat(2,2) = params[8];
    return mat;
}



