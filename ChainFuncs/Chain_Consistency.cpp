#include "../Chain.h"



//////////////////////////////////////////////////////////
////////////// CHAIN CONSISTENCY METHODS /////////////////
//////////////////////////////////////////////////////////

bool Chain::config_consistent() {
/*
    Check for consistency between triads and positions.
*/
    // check tangent position consistency
    for (unsigned i=0;i<num_bp-1;i++){
        if (arma::norm(bp_pos.col(i)+triads.slice(i).col(2)*disc_len - bp_pos.col(i+1)) >  MAX_POSITION_DISCREPANCY) {
            std::cout << "Configuration inconsistent " << i << " " << i+1 << ": " << arma::norm(bp_pos.col(i)+triads.slice(i).col(2)*disc_len - bp_pos.col(i+1)) << std::endl;
            return false;
        }
    }
    if (closed_topology) {
        if (arma::norm(bp_pos.col(num_bp-1)+triads.slice(num_bp-1).col(2)*disc_len - bp_pos.col(0)) >  MAX_POSITION_DISCREPANCY) {
            std::cout << "Configuration inconsistent " << num_bp-1 << " " << 0 << ": " <<  arma::norm(bp_pos.col(num_bp-1)+triads.slice(num_bp-1).col(2)*disc_len - bp_pos.col(0)) << std::endl;
            return false;
        }
    }

    double mddet = 0;
    double ddet;
    // check triads SO3 consistency
    for (unsigned i=0;i<num_bp;i++){
        ddet = std::abs(1-arma::det(triads.slice(i)));
        if (ddet > mddet) {
            mddet = ddet;
        }
        if (mddet > MAX_SO3_DETERMINANT_DISCREPANCY) {
            std::cout << "SO3 inconsistency: determinant differs from 1 by  " << mddet << std::endl;
            return false;
        }
    }
    return true;
}



void Chain::restore_consistency() {
/*
    Restores consistency between triads and positions, which might have been eroded due to
    an accumulation of numerical errors.
*/
    std::cout << "Restorning chain consistency" << std::endl;

    // Positions
    int first=0;
    if (num_bp%2==1) {
        first=1;
        arma::colvec tang = triads.slice(0).col(2);
        tang = tang/arma::norm(tang);
        bp_pos.col(1) = bp_pos.col(0 ) + tang*disc_len;
    }
    int pairs=num_bp/2;
    int id1,id2,id3;
    arma::colvec v,nv,w,o,no;
    double lv,lo,hlv;
    double discl_sq = disc_len*disc_len;
    for (int i=0;i<pairs;i++) {
        id1 = first+i*2;
        id2 = id1+1;
        id3 = pmod(id1+2,num_bp);
        v  = bp_pos.col(id3) - bp_pos.col(id1);
        w  = bp_pos.col(id2) - bp_pos.col(id1);
        lv = arma::norm(v);
        nv = v/lv;
        o  = w - arma::dot(nv,w)*nv;
        lo = arma::norm(o);
        if (lo < 1e-12) {
            o = triads.slice(id1).col(0);
            lo = arma::norm(o);
        }
        no = o/lo;
        hlv = lv/2;
        arma::colvec p2 = bp_pos.col(id2);
        bp_pos.col(id2) = bp_pos.col(id1) + v/2 + std::sqrt(discl_sq-hlv*hlv)*no;
    }
    // Triads
//    int ip1;
    arma::colvec e1,e2,e3;
    for (unsigned i=0;i<num_bp;i++) {
        e3 = bp_pos.col(pmod(i+1,num_bp)) - bp_pos.col(i);
        e3 = e3/arma::norm(e3);
        e1 = triads.slice(i).col(0);
        e1 = e1 - arma::dot(e1,e3)*e3;
        e1 = e1/arma::norm(e1);
        e2 = triads.slice(i).col(1);
        e2 = e2 - arma::dot(e2,e3)*e3;
        e2 = e2/arma::norm(e2);
        triads.slice(i).col(0) = e1;
        triads.slice(i).col(1) = e2;
        triads.slice(i).col(2) = e3;
    }
}

bool Chain::check_link_conservation() {
    double Wr, Tw, Lk, diff_Lk;
    Wr = cal_langowski_writhe_1a();
    Tw = cal_twist();
    Lk = Wr+Tw;
    diff_Lk = std::abs(Lk-dLK);
    if (diff_Lk > 1) {
        std::cout << "Linking Number not Conserved" << std::endl;
        std::cout << "Fixed Val:  " << dLK << std::endl;
        std::cout << "Calculated: " << Lk  << std::endl;
        return false;
    }
    return true;
}










