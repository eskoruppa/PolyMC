#ifndef __UNBOUND_INCLUDED__
#define __UNBOUND_INCLUDED__

#define UNBOUND_DEBUG

#define EQUAL_EPS 1e-10

/* --------------------------------------------*/
//  INCLUDE PAIRS
#include "Pair.h"
#include "PairStyles/Pair_LJ.h"


#define PAIR_TYPE_IDENFIER_LENNARD_JONES "LJ"

/* --------------------------------------------*/

#define TYPE_IDENTIFIER                 "type"
#define PAIR_IDENTIFIER                 "pair"
#define DEFAULT_TYPE_IDENTIFIER         "default_type"

#define INIT_TYPE_RANGE_SET_IDENFIER    "type_set_range"
#define INIT_TYPE_SET_IDENFIER          "type_set"


#define TYPE_CHANGE_MONOMER_ID 0
#define TYPE_CHANGE_CHANGE_TO  1


#include "AtomType.h"

#include "../Chain.h"
#include "../SO3Methods.h"
#include "../ExtraFuncs.h"

// include ExcludedVolume for definitions
#include "../ExcludedVolume/ExcludedVolume.h"

//// include input functions
#include "../Input/InputRead.h"

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


inline bool within_EV_unmoved(int movtype) {
    return within_EV_typeA(movtype);
}
inline bool within_EV_collective(int movtype) {
    return ( within_EV_typeB(movtype) || within_EV_typeD(movtype));
}
inline bool within_EV_individual(int movtype) {
    return ( within_EV_typeD(movtype) || within_EV_typeE(movtype));
}


class Unbound;

class Unbound {

protected:

    Chain *     chain;
//    arma::mat*  new_pos;
//    arma::cube* new_trds;
//    arma::mat*  old_pos;
//    arma::cube* old_trds;
//    int     num_bp;
//    double  disc_len;


    std::vector<TypeGroup*>  type_groups;
    std::vector<Pair*>       pairs;

    unsigned num_types;
    unsigned num_pairs;

    int              default_type = 0;
    std::vector<int> full_type_map;

    std::vector<int> backup_full_type_map;


    std::vector<std::vector<int>>   type_changes;
    bool                            type_change_acceptance_pending = false;


public:

    Unbound(Chain * ch);
    ~Unbound();

public:
/*
    These backups are for situations where the configuration has become inconsistent and
    the simulation reverts back to a previous, consistent marker.
*/
    void set_backup();
    void revert_to_backup();

public:
    void    submit_type_changes(const std::vector<std::vector<int>> & changes);
    double  eval_type_changes(const arma::mat * pos);
    void    set_type_changes(bool accept);

protected:
    void    assign_trial_elements();


public:
    double Eval_Delta_Energy(   const arma::mat * pos_new,
                                const arma::mat * pos_old,
                                const std::vector<arma::ivec> & moved
                                );

protected:
    int gen_type_groups_moved_intervals(const std::vector<arma::ivec> & moved);


public:
    int atom_type_exists(int type_id);
    int get_index_in_type_groups(int type_id);

protected:
    void assign_type_group_elements();


/* ------------------------------------------------*/
/*
    INIT
*/

public:
    void read_inputfile(const std::string & pairinteraction_file);

    void init_monomer_type_assignments_from_file(const std::string & typeset_file);
    void init_monomer_type_assignments(const std::vector<std::vector<int>> & range_sets,
                                                const std::vector<std::vector<int>> & individual_sets);

protected:
    void init_atoms(const std::vector<std::vector<std::string>> &);
    void init_pairs(const std::vector<std::vector<std::string>> &);

};





















#endif
