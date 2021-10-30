#ifndef __DUMP_FULL_COUP_MATRIX__
#define __DUMP_FULL_COUP_MATRIX__
#include "Dump.h"

/*
    TODO: Write to binary file rather than textfile.
*/


class Dump_Full_Coup_Matrix;

class Dump_Full_Coup_Matrix: public Dump
{

public:
    Dump_Full_Coup_Matrix(Chain * ch, const std::string& filename);
    ~Dump_Full_Coup_Matrix();

//    void prod_dump();
//    void final_dump();

    arma::mat params2mat(std::vector<double> params);
};




#endif
