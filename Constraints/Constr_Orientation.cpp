#include "Constr_Orientation.h"

Constr_Orientation::Constr_Orientation(Chain * ch,const std::vector<int> & constr_group) :
Constraint(ch,constr_group)
{
    init_fix();
}

Constr_Orientation::Constr_Orientation(Chain * ch,const std::vector<int> & constr_group, arma::mat * backup_pos, arma::cube * backup_triads) :
Constraint(ch,constr_group,backup_pos,backup_triads)
{
    init_fix();
}

Constr_Orientation::~Constr_Orientation() {}

void Constr_Orientation::init_fix() {
    name = "Constr_Orientation";

    constr_group = filter_double_entries(constr_group);

    arma::colvec euler;
    arma::mat    Ti;
    arma::mat    Tj;
    for (int i=0;i<constr_group.size()-1;i++) {
        for (int j=i+1;j<constr_group.size();j++) {
            Ti = triads->slice(constr_group[i]);
            Tj = triads->slice(constr_group[j]);
            euler = ExtractTheta( Ti.t() * Tj );
            pair_eulers.push_back(euler);
        }
    }
}


bool Constr_Orientation::check_constraint(const std::vector<arma::ivec>* moved) {
    arma::colvec euler;
    arma::mat    Ti;
    arma::mat    Tj;
    int id = -1;
    for (int i=0;i<constr_group.size()-1;i++) {
        for (int j=i+1;j<constr_group.size();j++) {
            id++;
            Ti = triads->slice(constr_group[i]);
            Tj = triads->slice(constr_group[j]);
            euler = ExtractTheta( Ti.t() * Tj );
            if (arma::norm(euler-pair_eulers[id]) > DIST_THRESHOLD) {
                return false;
            }
        }
    }
    return true;
}
