#include "EV_Beads.h"


ExcludedVolumeBead::ExcludedVolumeBead(int EVBead_id, int refbead_id, int num_contained, double rad, arma::mat* pos) :
id(EVBead_id), id_refbead(refbead_id), radius(rad), num_contained_bp(num_contained), all_pos(pos) {
    assert(num_contained_bp%2==0);
    int range = (num_contained_bp-1)/2;
    assert(refbead_id-range>=0);
    assert(refbead_id+range<all_pos->n_cols);
    for (int i=refbead_id-range;i<=refbead_id+range;i++) {
        contained_bp.push_back(i);
    }
}

void ExcludedVolumeBead::set_ProximityCluster(ProximityCluster* pc) {
    PC = pc;
}

int ExcludedVolumeBead::get_id() {
    return id;
}

int ExcludedVolumeBead::get_id_refbead(){
    return id_refbead;
}

int ExcludedVolumeBead::get_num_contained_bp() {
    return num_contained_bp;
}

std::vector<int>* ExcludedVolumeBead::get_contained_bp() {
    return &contained_bp;
}





