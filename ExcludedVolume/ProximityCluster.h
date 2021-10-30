#ifndef __PROXCLUSTER__
#define __PROXCLUSTER__

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

class ProximityCluster;

class ProximityCluster {

protected:
    double max_displacement;
    double radius;
    double interaction_distance;

//    vector<ExcludedVolumeBead*> contained;
    int ref_bp;
    int first;
    int last;
    int num_contained;

    int ID;
    int num_PC;
    std::vector<ProximityCluster*> PC_list;
    std::vector<ProximityCluster*> proximity_PCptrs;
    std::vector<int>          proximity_PCids;

    std::vector<arma::ivec> proximity_interval;
    int                proximity_interval_self;
    int                num_proximity_intervals;

    arma::colvec pos;
    arma::colvec prev_pos;

public:

    ProximityCluster(int EV_first, int EV_last, int reference_bp, const arma::colvec& init_pos, double radius);
    void set_list(std::vector<ProximityCluster*> bc_list, int own_ID);
    int  get_ID();

    int get_ref_bp();
    int get_first();
    int get_last();

    bool          check_displacement(arma::colvec& new_pos);
    void          set_pos(arma::colvec& new_pos);
    arma::colvec* get_pos();

    void clear_proximity_list();
    void calculate_proximity_ascending();
    void calculate_proximity_full();
protected:
    void add_to_proximity(int id);
public:
    void build_proximity_interval();
    int  get_proximity_interval( std::vector<arma::ivec>** proximity_interval_ptr);

    void get_proximity_tail(std::vector<arma::ivec>** proximity_interval_ptr,
                                arma::colvec& rel_interval,
                                int left_limit,
                                int right_limit);



};

#endif
