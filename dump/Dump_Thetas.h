#ifndef __DUMP_THETAS_H__
#define __DUMP_THETAS_H__
#include "Dump.h"


////////////////////////////////////////////////////////////////////////////////
/*
    Dump all Thetas every N_dump steps
*/

class Dump_Thetas;

class Dump_Thetas: public Dump
{
protected:
    int counter;

public:
    Dump_Thetas(Chain * ch, int N_dump, const std::string& filename,bool append=false);
    ~Dump_Thetas();

    void prod_dump();
};


////////////////////////////////////////////////////////////////////////////////
/*
    Dump average of all Thetas, calculated every N_dump steps
*/

class Dump_avgThetas;

class Dump_avgThetas: public Dump
{
protected:
    int counter;
    arma::mat       Thetas;
    arma::mat       Thetas_sq;
    arma::colvec    kappa;

public:
    Dump_avgThetas(Chain * ch, int N_dump, const std::string& filename,bool append=false);
    ~Dump_avgThetas();

    void prod_dump();
    void final_dump();
};


#endif
