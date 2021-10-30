#ifndef __INCLUDED_CONSTR_ORIENTATION__
#define __INCLUDED_CONSTR_ORIENTATION__

#include "Constraint.h"

#define DIST_THRESHOLD 1e-10

class Constr_Orientation;

class Constr_Orientation: public Constraint {
protected:
    std::vector<arma::colvec> pair_eulers;

public:

    Constr_Orientation(Chain * ch,const std::vector<int> & constr_group);
    Constr_Orientation(Chain * ch,const std::vector<int> & constr_group, arma::mat * backup_pos, arma::cube * backup_triads);
    ~Constr_Orientation();

protected:
    void init_fix();


public:
    bool check_constraint(const std::vector<arma::ivec>* moved);

};

#endif
