#ifndef __INCLUDED_CONSTRAINT__
#define __INCLUDED_CONSTRAINT__

#include "../Chain.h"
#include "../SO3Methods.h"
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


class Constraint;

class Constraint {
protected:
    std::string name = "Constraint";

    std::vector<int> constr_group;

    Chain * chain;
    arma::mat*  bp_pos;
    arma::cube* triads;

    int num_bp;
    int num_bps;
    double disc_len;

    int counter;
    int counter_reject;

    arma::mat  * ptr_bp_pos_backup;
    arma::cube * ptr_triads_backup;

    arma::mat  bp_pos_backup;
    arma::cube triads_backup;

    bool intrinsic_backup = false;
    bool set_backup = false;


public:

    Constraint(Chain * ch,const std::vector<int> & constr_group);
    Constraint(Chain * ch,const std::vector<int> & constr_group, arma::mat * backup_pos, arma::cube * backup_triads);
    virtual ~Constraint();

    bool check_valid_constr_group();

    arma::mat*  get_bp_pos_backup();
    arma::cube* get_triads_backup();
    void        set_backup_conf(arma::mat * backup_pos, arma::cube * backup_triads);

public:
    bool check(const std::vector<arma::ivec>* moved);
    virtual bool check_constraint(const std::vector<arma::ivec>* moved);


protected:
    void    transfer_config(const std::vector<arma::ivec>* moved, arma::mat* bp_pos_from, arma::mat* bp_pos_to, arma::cube* triads_from, arma::cube* triads_to );
    const std::vector<int> filter_double_entries(const std::vector<int> & group);

public:
    void    revert_to_backup();
    void    set_current_as_backup();

public:
    double rejection_rate();

// setters
public:
    void allow_setting_backup(bool allowed=true);


};

#endif
