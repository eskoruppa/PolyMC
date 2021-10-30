#include "Constr_Position.h"

Constr_Position::Constr_Position(Chain * ch, const std::vector<int> & constr_group) :
Constraint(ch,constr_group)
{
    init_fix();
}

Constr_Position::Constr_Position(Chain * ch,const std::vector<int> & constr_group, arma::mat * backup_pos, arma::cube * backup_triads) :
Constraint(ch,constr_group,backup_pos,backup_triads)
{
    init_fix();
}

Constr_Position::~Constr_Position() {}

void Constr_Position::init_fix() {
    name = "Constr_Position";

    constr_group = filter_double_entries(constr_group);
    double dist;
    for (int i=0;i<constr_group.size()-1;i++) {
//        int j = constr_group.size()-1;
        for (int j=i+1;j<constr_group.size();j++) {
            dist = arma::norm(bp_pos->col(constr_group[i])-bp_pos->col(constr_group[j]));
            pair_dists.push_back(dist);
        }
    }
}


bool Constr_Position::check_constraint(const std::vector<arma::ivec>* moved) {
    double dist;
    int id = -1;
    for (int i=0;i<constr_group.size()-1;i++) {
//        int j = constr_group.size()-1;
        for (int j=i+1;j<constr_group.size();j++) {
            id++;
            dist = arma::norm(bp_pos->col(constr_group[i])-bp_pos->col(constr_group[j]));
            if (std::abs(dist-pair_dists[id]) > DIST_THRESHOLD) {
                return false;
            }
        }
    }
//    std::cout << constr_group.size() << std::endl;

    return true;
}
