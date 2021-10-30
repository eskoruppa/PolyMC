#ifndef __INCLUDED_CONSTR_POS__
#define __INCLUDED_CONSTR_POS__

#include "Constraint.h"

#define DIST_THRESHOLD 1e-10

class Constr_Position;

class Constr_Position: public Constraint {
protected:
    std::vector<double> pair_dists;

public:

    Constr_Position(Chain * ch,const std::vector<int> & constr_group);
    Constr_Position(Chain * ch,const std::vector<int> & constr_group, arma::mat * backup_pos, arma::cube * backup_triads);
    ~Constr_Position();

protected:
    void init_fix();


public:
    bool check_constraint(const std::vector<arma::ivec>* moved);

};

#endif
