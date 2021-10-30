#ifndef __DUMP_H__
#define __DUMP_H__

#include "../Chain.h"
#include "../SO3Methods.h"

#include <armadillo> // armadillo
#include <iostream>  // input output
#include <iomanip>   // std::setprecision
#include <string>    // string
#include <algorithm>
#include <stdexcept> // exceptions
#include <fstream>   // read file
#include <cmath>     // math lib
#include <random>    // random lib


class Dump;

class Dump {

protected:
    Chain * chain;
    std::vector<BPStep*> BPS;

    arma::cube * triads;
    arma::mat  * pos;

    int num_bp;     // number base pairs
    int num_bps;    // number base pair steps
    std::string fn; // filename
    int  N_step;    // dump at every Nth step
    int  step;      // step count
    bool app;       // append
    double disc_len;

public:

    Dump(Chain * ch, int N_dump, const std::string& filename,bool append=false);
    virtual ~Dump();

    void dump();
    std::string get_filename();

//    void init();
    virtual void prod_dump();
    virtual void final_dump();

protected:
    bool fexists(const std::string& filename);

};

#endif
