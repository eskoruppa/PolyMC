#ifndef __PAIR_INCLUDED__
#define __PAIR_INCLUDED__

#include "../Chain.h"
#include "AtomType.h"
#include "../ExtraFuncs.h"

#define _USE_MATH_DEFINES

#include <armadillo> // armadillo
#include <iostream>  // input output
#include <string>    // string
#include <algorithm>
#include <stdexcept> // exceptions
#include <fstream>   // read file
#include <cmath>     // math lib
#include <vector>    // vector
#include <cassert>   // assert


#define PAIR_INFINITY 1e10


#define PAIR_DEBUG



#define TEST \
    for (int i=0;i<10;i++) { \
        std::cout << i << std::endl; \
    }





class Pair;

class Pair {

protected:

    Chain * chain;
    int     num_bp;
    double  disc_len;

    std::vector<double> params;

    bool single_group;

    TypeGroup * type_group_A;
    TypeGroup * type_group_B;

    AtomType * type_A;
    AtomType * type_B;

    int type_id_A;
    int type_id_B;

    std::vector<int> * group_A;
    std::vector<int> * group_B;

    double lower_cutoff; // defines hard repulsion zone
    double upper_cutoff; // cutoff for the potential

    std::vector<int> * A_unmoved;
    std::vector<int> * A_moved;
    std::vector<int> * A_individual;
    std::vector<int> * B_unmoved;
    std::vector<int> * B_moved;
    std::vector<int> * B_individual;

    int * num_moved_A;
    int * num_moved_B;

//    arma::mat*  new_pos;
//    arma::cube* new_trds;
//    arma::mat*  old_pos;
//    arma::cube* old_trds;

public:

    Pair(   Chain * ch,
            TypeGroup * type_group_A,
            TypeGroup * type_group_B,
            double cutoff,
            std::vector<double> params);

    Pair();
    virtual ~Pair();

    /*---------------------------------------------------------------------------------------------------------------------*/

    virtual double Eval_change_type(const arma::mat * pos);
    virtual double Eval_change_type(const arma::mat * pos,
                            int monomer_id,
                            int type_from,
                            int type_to);

    /*---------------------------------------------------------------------------------------------------------------------*/

    virtual double Eval_Energy(const arma::mat * pos);

    virtual double Eval_Delta_Energy(   const arma::mat * pos_new,
                                        const arma::mat * pos_old);
    virtual double Eval_Delta_Energy(   const arma::mat * pos_new,
                                        const arma::mat * pos_old,
                                        const std::vector<int> * A_unmoved,
                                        const std::vector<int> * A_moved,
                                        const std::vector<int> * A_individual,
                                        const std::vector<int> * B_unmoved,
                                        const std::vector<int> * B_moved,
                                        const std::vector<int> * B_individual);

    /*---------------------------------------------------------------------------------------------------------------------*/

protected:

    double Potential(double distance);

    unsigned find_id_in_group_A(int monomer_id);
    unsigned find_id_in_group_B(int monomer_id);

};





















#endif
