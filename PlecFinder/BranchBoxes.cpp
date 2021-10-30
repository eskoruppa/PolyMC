#include "BranchBoxes.h"

BranchBox::BranchBox(const std::vector<int> & xlim,const std::vector<int> & ylim, arma::mat * WM) {
    set_lims(xlim,ylim);
    this->WM   = WM;
}

void    BranchBox::set_WM(arma::mat * WM) {
    this->WM = WM;
}

void    BranchBox::set_xlim(const std::vector<int> & xlim) {
    this->xlim   = xlim;
    check_finite_size();
}

void    BranchBox::set_ylim(const std::vector<int> & ylim) {
    this->ylim   = ylim;
    check_finite_size();
}

void    BranchBox::set_lims(const std::vector<int> & xlim,const std::vector<int> & ylim) {
    this->xlim   = xlim;
    this->ylim   = ylim;
    check_finite_size();
}

double  BranchBox::get_writhe() {
    if (!active) {
        writhe = 0;
        writhe_calculated=false;
        return 0;
    }
    if (!writhe_calculated) {
        writhe = arma::accu((*WM)(arma::span(xlim[0],xlim[1]-1),arma::span(ylim[0],ylim[1]-1)))*2;
        writhe_calculated = true;
    }
    if (writhe <= 0) active = false;
    return writhe;
}

bool    BranchBox::is_active() {
    return active;
}

void BranchBox::check_finite_size() {
    if (xlim[0] >= xlim[1] || ylim[0] >= ylim[1] )
        active = false;
    else {
        active = true;
    }
    writhe_calculated = false;
}
