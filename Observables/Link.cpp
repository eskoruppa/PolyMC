#include "../Chain.h"

//changed

double Chain::cal_quick_writhe() {
    if (link_constrained) {
        double Twist = cal_twist(0,num_bp);
        return dLK - Twist;
    }
    else {
        return cal_langowski_writhe_1a(0.05);
    }
}

double Chain::cal_fuller_writhe( int from, int to, arma::colvec& z_dir) {
    if (from <= 0) from = 1;
    if (to >= num_bp) to = num_bp-1;
    double Wr = 0;
    double denom;
    arma::colvec e3p1,e3m1,e3;
    for (int i=from;i<to;i++) {
        e3m1 = triads.slice(i-1).col(2);
        e3p1 = triads.slice(i+1).col(2);
        e3   = triads.slice( i ).col(2);
        denom = 1+dot(e3,z_dir);
        Wr = Wr + dot( z_dir, cross ( e3,(e3p1-e3m1)*0.5) ) / denom ;
    }
    return Wr/(2*M_PI);
}

double Chain::cal_twist(int from, int to) {
    if (from < 0) from = 0;
    if (to >= num_bps) to = num_bps-1;
    double Tw = 0;
    arma::colvec Theta;
    for (int i=from;i<to;i++) {
        Theta = *BPS[i]->get_Theta();
        Tw = Tw + Theta(2);
    }
    return Tw/(2*M_PI);
}

double Chain::cal_twist() {
    return cal_twist(0,num_bps);
}

double Chain::gauss_writhe(double density) {

    double Wr = 0;
    int          id1,id1p1,id1m1,id2,id2p1,id2m1;
    arma::colvec r1,r2,dr1,dr2,r12;
    double       len_r12;

    int num_points = (int)(num_bp*density);
    int point_dist = num_bp/num_points;

    double da = point_dist;

    for (int i=0;i<num_points-1;i++) {
        for (int j=i+1;j<num_points;j++) {
            id1   = i*point_dist;
            id1p1 = (i+1)*point_dist;
            id1m1 = pmod((i-1)*point_dist,num_bp);
            id2   = j*point_dist;
            id2p1 = pmod((j+1)*point_dist,num_bp);
            id2m1 = (j-1)*point_dist;

            r1 = bp_pos.col(id1);
            r2 = bp_pos.col(id2);
//            dr1 = (bp_pos.col(id1p1)-bp_pos.col(id1m1))/da;
//            dr2 = (bp_pos.col(id2p1)-bp_pos.col(id2m1))/da;
            dr1 = (bp_pos.col(id1p1)-r1)/da;
            dr2 = (bp_pos.col(id2p1)-r2)/da;
            r12 = r2-r1;
            len_r12 = arma::norm(r12);

            Wr += arma::dot(arma::cross(dr1,dr2),r12)/(len_r12*len_r12*len_r12);
        }
    }
    Wr /= -2*M_PI;
    return Wr;
}


double Chain::cal_langowski_writhe_1a() {
    double Wr   = 0;
    int    id1,id2,id3,id4;
    arma::colvec r12,r13,r14,r23,r24,r34;
    arma::colvec n1,n2,n3,n4;
    double Omega;
    double nn1,nn2,nn3,nn4;

    for (unsigned i=1;i<num_bp-1;i++) {
        id3 = i;
        id4 = i+1;
        r34 = bp_pos.col(id4) - bp_pos.col(id3);
        for (int j=0;j<i-1;j++) {
            id1 = j;
            id2 = j+1;
            r12 = bp_pos.col(id2) - bp_pos.col(id1);
            r13 = bp_pos.col(id3) - bp_pos.col(id1);
            r14 = bp_pos.col(id4) - bp_pos.col(id1);
            r23 = bp_pos.col(id3) - bp_pos.col(id2);
            r24 = bp_pos.col(id4) - bp_pos.col(id2);

            n1  = arma::cross(r13,r14);
            nn1 = arma::norm(n1);
            if (nn1 > 1e-10) {
                n1 = n1/nn1;
//                continue;
            }
            n2  = arma::cross(r14,r24);
            nn2 = arma::norm(n2);
            if (nn2 > 1e-10) {
                n2 = n2/nn2;
//                continue;
            }
            n3  = arma::cross(r24,r23);
            nn3 = arma::norm(n3);
            if (nn3 > 1e-10) {
                n3 = n3/nn3;
//                continue;
            }
            n4  = arma::cross(r23,r13);
            nn4 = arma::norm(n4);
            if (nn4 > 1e-10) {
                n4 = n4/nn4;
//                continue;
            }

            Omega =  std::asin(arma::dot(n1,n2));
            Omega += std::asin(arma::dot(n2,n3));
            Omega += std::asin(arma::dot(n3,n4));
            Omega += std::asin(arma::dot(n4,n1));
            Omega = Omega * sgn(arma::dot(cross(r34,r12),r13));
            if (Omega == Omega) {
                Wr = Wr + Omega;
            }

        }
    }
    if (closed_topology == true) {
        id3 = num_bp-1;
        id4 = 0;
        r34 = bp_pos.col(id4) - bp_pos.col(id3);
        for (int j=1;j<id3-1;j++) {
            id1 = j;
            id2 = j+1;
            r12 = bp_pos.col(id2) - bp_pos.col(id1);
            r13 = bp_pos.col(id3) - bp_pos.col(id1);
            r14 = bp_pos.col(id4) - bp_pos.col(id1);
            r23 = bp_pos.col(id3) - bp_pos.col(id2);
            r24 = bp_pos.col(id4) - bp_pos.col(id2);

            n1  = arma::cross(r13,r14);
            nn1 = arma::norm(n1);
            if (nn1 > 1e-10) {
                n1 = n1/nn1;
//                continue;
            }
            n2  = arma::cross(r14,r24);
            nn2 = arma::norm(n2);
            if (nn2 > 1e-10) {
                n2 = n2/nn2;
//                continue;
            }
            n3  = arma::cross(r24,r23);
            nn3 = arma::norm(n3);
            if (nn3 > 1e-10) {
                n3 = n3/nn3;
//                continue;
            }
            n4  = arma::cross(r23,r13);
            nn4 = arma::norm(n4);
            if (nn4 > 1e-10) {
                n4 = n4/nn4;
//                continue;
            }

            Omega =  std::asin(arma::dot(n1,n2));
            Omega += std::asin(arma::dot(n2,n3));
            Omega += std::asin(arma::dot(n3,n4));
            Omega += std::asin(arma::dot(n4,n1));
            Omega = Omega * sgn(arma::dot(cross(r34,r12),r13));
            if (Omega == Omega) {
                Wr = Wr + Omega;
            }

        }
    }
    return Wr/(2*M_PI);
}

double Chain::cal_langowski_writhe_1a(double density) {

    double Wr   = 0;
    int    id1,id2,id3,id4;
    arma::colvec r12,r13,r14,r23,r24,r34;
    arma::colvec n1,n2,n3,n4;
    double Omega;
    double nn1,nn2,nn3,nn4;

    int num_points = (int)(num_bp*density);
    int point_dist = num_bp/num_points;

    for (int i=1;i<num_points-1;i++) {
        id3 = i*point_dist;
        id4 = (i+1)*point_dist;
        r34 = bp_pos.col(id4) - bp_pos.col(id3);
        for (int j=0;j<i-1;j++) {
            id1 = j*point_dist;
            id2 = (j+1)*point_dist;
            r12 = bp_pos.col(id2) - bp_pos.col(id1);
            r13 = bp_pos.col(id3) - bp_pos.col(id1);
            r14 = bp_pos.col(id4) - bp_pos.col(id1);
            r23 = bp_pos.col(id3) - bp_pos.col(id2);
            r24 = bp_pos.col(id4) - bp_pos.col(id2);

            n1  = arma::cross(r13,r14);
            nn1 = arma::norm(n1);
            if (nn1 > 1e-10) {
                n1 = n1/nn1;
//                continue;
            }
            n2  = arma::cross(r14,r24);
            nn2 = arma::norm(n2);
            if (nn2 > 1e-10) {
                n2 = n2/nn2;
//                continue;
            }
            n3  = arma::cross(r24,r23);
            nn3 = arma::norm(n3);
            if (nn3 > 1e-10) {
                n3 = n3/nn3;
//                continue;
            }
            n4  = arma::cross(r23,r13);
            nn4 = arma::norm(n4);
            if (nn4 > 1e-10) {
                n4 = n4/nn4;
//                continue;
            }

            Omega =  std::asin(arma::dot(n1,n2));
            Omega += std::asin(arma::dot(n2,n3));
            Omega += std::asin(arma::dot(n3,n4));
            Omega += std::asin(arma::dot(n4,n1));
            Omega = Omega * sgn(arma::dot(cross(r34,r12),r13));
            if (Omega == Omega) {
                Wr = Wr + Omega;
            }

        }
    }
    if (closed_topology == true) {
        id3 = (num_points-1)*point_dist;
        id4 = 0;
        r34 = bp_pos.col(id4) - bp_pos.col(id3);
        for (int j=1;j<(num_points-1)-1;j++) {
            id1 = j*point_dist;
            id2 = (j+1)*point_dist;
            r12 = bp_pos.col(id2) - bp_pos.col(id1);
            r13 = bp_pos.col(id3) - bp_pos.col(id1);
            r14 = bp_pos.col(id4) - bp_pos.col(id1);
            r23 = bp_pos.col(id3) - bp_pos.col(id2);
            r24 = bp_pos.col(id4) - bp_pos.col(id2);

            n1  = arma::cross(r13,r14);
            nn1 = arma::norm(n1);
            if (nn1 > 1e-10) {
                n1 = n1/nn1;
//                continue;
            }
            n2  = arma::cross(r14,r24);
            nn2 = arma::norm(n2);
            if (nn2 > 1e-10) {
                n2 = n2/nn2;
//                continue;
            }
            n3  = arma::cross(r24,r23);
            nn3 = arma::norm(n3);
            if (nn3 > 1e-10) {
                n3 = n3/nn3;
//                continue;
            }
            n4  = arma::cross(r23,r13);
            nn4 = arma::norm(n4);
            if (nn4 > 1e-10) {
                n4 = n4/nn4;
//                continue;
            }

            Omega =  std::asin(arma::dot(n1,n2));
            Omega += std::asin(arma::dot(n2,n3));
            Omega += std::asin(arma::dot(n3,n4));
            Omega += std::asin(arma::dot(n4,n1));
            Omega = Omega * sgn(arma::dot(cross(r34,r12),r13));
            if (Omega == Omega) {
                Wr = Wr + Omega;
            }

        }
    }
    return Wr/(2*M_PI);
}

double Chain::cal_langowski_writhe_1a(const arma::mat& pos, bool closed) {
    int segs = pos.n_cols;
    double Wr   = 0;

    int    id1,id2,id3,id4;
    arma::colvec r12,r13,r14,r23,r24,r34;
    arma::colvec n1,n2,n3,n4;
    double Omega;
    double nn1,nn2,nn3,nn4;

    for (int i=1;i<segs-1;i++) {

        id3 = i;
        id4 = i+1;
        r34 = pos.col(id4) - pos.col(id3);
        for (int j=0;j<i-1;j++) {
            id1 = j;
            id2 = j+1;
            r12 = pos.col(id2) - pos.col(id1);
            r13 = pos.col(id3) - pos.col(id1);
            r14 = pos.col(id4) - pos.col(id1);
            r23 = pos.col(id3) - pos.col(id2);
            r24 = pos.col(id4) - pos.col(id2);

            n1  = arma::cross(r13,r14);
            nn1 = arma::norm(n1);
            if (nn1 > 1e-10) {
                n1 = n1/nn1;
//                continue;
            }
            n2  = arma::cross(r14,r24);
            nn2 = arma::norm(n2);
            if (nn2 > 1e-10) {
                n2 = n2/nn2;
//                continue;
            }
            n3  = arma::cross(r24,r23);
            nn3 = arma::norm(n3);
            if (nn3 > 1e-10) {
                n3 = n3/nn3;
//                continue;
            }
            n4  = arma::cross(r23,r13);
            nn4 = arma::norm(n4);
            if (nn4 > 1e-10) {
                n4 = n4/nn4;
//                continue;
            }

            Omega =  std::asin(arma::dot(n1,n2));
            Omega += std::asin(arma::dot(n2,n3));
            Omega += std::asin(arma::dot(n3,n4));
            Omega += std::asin(arma::dot(n4,n1));
            Omega = Omega * sgn(arma::dot(cross(r34,r12),r13));
            if (Omega == Omega) {
                Wr = Wr + Omega;
            }
        }
    }
    if (closed == true) {
        id3 = segs-1;
        id4 = 0;
        r34 = pos.col(id4) - pos.col(id3);
        for (int j=1;j<id3-1;j++) {
            id1 = j;
            id2 = j+1;
            r12 = pos.col(id2) - pos.col(id1);
            r13 = pos.col(id3) - pos.col(id1);
            r14 = pos.col(id4) - pos.col(id1);
            r23 = pos.col(id3) - pos.col(id2);
            r24 = pos.col(id4) - pos.col(id2);

            n1  = arma::cross(r13,r14);
            nn1 = arma::norm(n1);
            if (nn1 > 1e-10) {
                n1 = n1/nn1;
//                continue;
            }
            n2  = arma::cross(r14,r24);
            nn2 = arma::norm(n2);
            if (nn2 > 1e-10) {
                n2 = n2/nn2;
//                continue;
            }
            n3  = arma::cross(r24,r23);
            nn3 = arma::norm(n3);
            if (nn3 > 1e-10) {
                n3 = n3/nn3;
//                continue;
            }
            n4  = arma::cross(r23,r13);
            nn4 = arma::norm(n4);
            if (nn4 > 1e-10) {
                n4 = n4/nn4;
//                continue;
            }

            Omega =  std::asin(arma::dot(n1,n2));
            Omega += std::asin(arma::dot(n2,n3));
            Omega += std::asin(arma::dot(n3,n4));
            Omega += std::asin(arma::dot(n4,n1));
            Omega = Omega * sgn(arma::dot(cross(r34,r12),r13));
            if (Omega == Omega) {
                Wr = Wr + Omega;
            }

        }
    }
    return Wr/(2*M_PI);
}

void Chain::langowski_writhe_elements(arma::mat* writhe_elem) {
    *writhe_elem = arma::zeros(num_bp,num_bp);

    arma::colvec r12,r13,r14,r23,r24,r34;
    arma::colvec n1,n2,n3,n4;
    double Omega;
    double fac4pi = 1./(4*M_PI);
    double vecnorm;

    /*
        Calculate Gauss integral segment contributions except those involving the last segment
    */

    for (unsigned i=0;i<num_bp-3;i++) {
        unsigned ip1 = i+1;
        r12 = bp_pos.col(ip1)-bp_pos.col(i);

        for (unsigned j=i+2;j<num_bp-1;j++) {
            unsigned jp1 = j+1;

            r13 = bp_pos.col(j)   - bp_pos.col(i);
            r14 = bp_pos.col(jp1) - bp_pos.col(i);
            r23 = bp_pos.col(j)   - bp_pos.col(ip1);
            r24 = bp_pos.col(jp1) - bp_pos.col(ip1);
            r34 = bp_pos.col(jp1) - bp_pos.col(j);

            n1      = arma::cross(r13,r14);
            vecnorm = arma::norm(n1);
            if (vecnorm > 1e-10) {
                n1 = n1/vecnorm;
            }

            n2      = arma::cross(r14,r24);
            vecnorm = arma::norm(n2);
            if (vecnorm > 1e-10) {
                n2 = n2/vecnorm;
            }

            n3      = arma::cross(r24,r23);
            vecnorm = arma::norm(n3);
            if (vecnorm > 1e-10) {
                n3 = n3/vecnorm;
            }

            n4      = arma::cross(r23,r13);
            vecnorm = arma::norm(n4);
            if (vecnorm > 1e-10) {
                n4 = n4/vecnorm;
            }

            Omega =  std::asin(arma::dot(n1,n2));
            Omega += std::asin(arma::dot(n2,n3));
            Omega += std::asin(arma::dot(n3,n4));
            Omega += std::asin(arma::dot(n4,n1));
            Omega = Omega*fac4pi*sgn( arma::dot( arma::cross(r34,r12),r13) );

            if (Omega==Omega) {
                (*writhe_elem)(i,j) = Omega;
                (*writhe_elem)(j,i) = Omega;
            }
        }
    }

    /*
        Calculate Gauss integral segment contributions involving the last segment
    */

    unsigned j = num_bp-1;
    arma::colvec lastpos = bp_pos.col(j)+triads.slice(j).col(2)*disc_len;
    for (unsigned i=1;i<num_bp-2;i++) {
        unsigned ip1 = i+1;

        r12 = bp_pos.col(ip1) - bp_pos.col(i);
        r13 = bp_pos.col(j)   - bp_pos.col(i);
        r14 = lastpos         - bp_pos.col(i);
        r23 = bp_pos.col(j)   - bp_pos.col(ip1);
        r24 = lastpos         - bp_pos.col(ip1);
        r34 = lastpos         - bp_pos.col(j);

        n1      = arma::cross(r13,r14);
        vecnorm = arma::norm(n1);
        if (vecnorm > 1e-10) {
            n1 = n1/vecnorm;
        }

        n2      = arma::cross(r14,r24);
        vecnorm = arma::norm(n2);
        if (vecnorm > 1e-10) {
            n2 = n2/vecnorm;
        }

        n3      = arma::cross(r24,r23);
        vecnorm = arma::norm(n3);
        if (vecnorm > 1e-10) {
            n3 = n3/vecnorm;
        }

        n4      = arma::cross(r23,r13);
        vecnorm = arma::norm(n4);
        if (vecnorm > 1e-10) {
            n4 = n4/vecnorm;
        }

        Omega =  std::asin(arma::dot(n1,n2));
        Omega += std::asin(arma::dot(n2,n3));
        Omega += std::asin(arma::dot(n3,n4));
        Omega += std::asin(arma::dot(n4,n1));
        Omega = Omega*fac4pi*sgn( arma::dot( arma::cross(r34,r12),r13) );

        if (Omega==Omega) {
            (*writhe_elem)(i,j) = Omega;
            (*writhe_elem)(j,i) = Omega;
        }
    }
}

void Chain::langowski_writhe_elements(arma::mat* writhe_elem, bool closed) {

    if (closed) *writhe_elem = arma::zeros(num_bp,num_bp);
    else        *writhe_elem = arma::zeros(num_bp-1,num_bp-1);

    arma::colvec r12,r13,r14,r23,r24,r34;
    arma::colvec n1,n2,n3,n4;
    double Omega;
    double fac4pi = 1./(4*M_PI);
    double vecnorm;

    /*
        Calculate Gauss integral segment contributions except those involving the last segment
    */

    for (unsigned i=0;i<num_bp-3;i++) {
        unsigned ip1 = i+1;
        r12 = bp_pos.col(ip1)-bp_pos.col(i);

        for (unsigned j=i+2;j<num_bp-1;j++) {
            unsigned jp1 = j+1;

            r13 = bp_pos.col(j)   - bp_pos.col(i);
            r14 = bp_pos.col(jp1) - bp_pos.col(i);
            r23 = bp_pos.col(j)   - bp_pos.col(ip1);
            r24 = bp_pos.col(jp1) - bp_pos.col(ip1);
            r34 = bp_pos.col(jp1) - bp_pos.col(j);

            n1      = arma::cross(r13,r14);
            vecnorm = arma::norm(n1);
            if (vecnorm > 1e-10) {
                n1 = n1/vecnorm;
            }

            n2      = arma::cross(r14,r24);
            vecnorm = arma::norm(n2);
            if (vecnorm > 1e-10) {
                n2 = n2/vecnorm;
            }

            n3      = arma::cross(r24,r23);
            vecnorm = arma::norm(n3);
            if (vecnorm > 1e-10) {
                n3 = n3/vecnorm;
            }

            n4      = arma::cross(r23,r13);
            vecnorm = arma::norm(n4);
            if (vecnorm > 1e-10) {
                n4 = n4/vecnorm;
            }

            Omega =  std::asin(arma::dot(n1,n2));
            Omega += std::asin(arma::dot(n2,n3));
            Omega += std::asin(arma::dot(n3,n4));
            Omega += std::asin(arma::dot(n4,n1));
            Omega = Omega*fac4pi*sgn( arma::dot( arma::cross(r34,r12),r13) );

            if (Omega==Omega) {
                (*writhe_elem)(i,j) = Omega;
                (*writhe_elem)(j,i) = Omega;
            }
        }
    }

    /*
        Calculate Gauss integral segment contributions involving the last segment
    */
    if (closed) {
        unsigned j = num_bp-1;
        arma::colvec lastpos = bp_pos.col(j)+triads.slice(j).col(2)*disc_len;
        for (unsigned i=1;i<num_bp-2;i++) {
            unsigned ip1 = i+1;

            r12 = bp_pos.col(ip1) - bp_pos.col(i);
            r13 = bp_pos.col(j)   - bp_pos.col(i);
            r14 = lastpos         - bp_pos.col(i);
            r23 = bp_pos.col(j)   - bp_pos.col(ip1);
            r24 = lastpos         - bp_pos.col(ip1);
            r34 = lastpos         - bp_pos.col(j);

            n1      = arma::cross(r13,r14);
            vecnorm = arma::norm(n1);
            if (vecnorm > 1e-10) {
                n1 = n1/vecnorm;
            }

            n2      = arma::cross(r14,r24);
            vecnorm = arma::norm(n2);
            if (vecnorm > 1e-10) {
                n2 = n2/vecnorm;
            }

            n3      = arma::cross(r24,r23);
            vecnorm = arma::norm(n3);
            if (vecnorm > 1e-10) {
                n3 = n3/vecnorm;
            }

            n4      = arma::cross(r23,r13);
            vecnorm = arma::norm(n4);
            if (vecnorm > 1e-10) {
                n4 = n4/vecnorm;
            }

            Omega =  std::asin(arma::dot(n1,n2));
            Omega += std::asin(arma::dot(n2,n3));
            Omega += std::asin(arma::dot(n3,n4));
            Omega += std::asin(arma::dot(n4,n1));
            Omega = Omega*fac4pi*sgn( arma::dot( arma::cross(r34,r12),r13) );

            if (Omega==Omega) {
                (*writhe_elem)(i,j) = Omega;
                (*writhe_elem)(j,i) = Omega;
            }
        }
    }
}

void Chain::langowski_writhe_elements(arma::mat* writhe_elem, double seg_size) {
    /*
        Calculates the writhe for segments sizes as specified in the argument. If the specified size is
        smaller than the discreitzation length it will be set to the discretization length.
    */

    if (seg_size<disc_len) seg_size=disc_len;

    int num_points = (int)(num_bp*disc_len/seg_size);
    int point_dist = num_bp/num_points;

    *writhe_elem = arma::zeros(num_points,num_points);

    arma::colvec r12,r13,r14,r23,r24,r34;
    arma::colvec n1,n2,n3,n4;
    double Omega;
    double fac4pi = 1./(4*M_PI);
    double vecnorm;

    /*
        Calculate Gauss integral segment contributions except those involving the last segment
    */

    int i,j;
    int ip1,jp1;

    for (int k=0;k<num_points-3;k++) {
        i   = k*point_dist;
        ip1 = (k+1)*point_dist;
        r12 = bp_pos.col(ip1)-bp_pos.col(i);

        for (int l=k+2;l<num_points-1;l++) {
            j   = l*point_dist;
            jp1 = (l+1)*point_dist;

            r13 = bp_pos.col(j)   - bp_pos.col(i);
            r14 = bp_pos.col(jp1) - bp_pos.col(i);
            r23 = bp_pos.col(j)   - bp_pos.col(ip1);
            r24 = bp_pos.col(jp1) - bp_pos.col(ip1);
            r34 = bp_pos.col(jp1) - bp_pos.col(j);

            n1      = arma::cross(r13,r14);
            vecnorm = arma::norm(n1);
            if (vecnorm > 1e-10) {
                n1 = n1/vecnorm;
            }

            n2      = arma::cross(r14,r24);
            vecnorm = arma::norm(n2);
            if (vecnorm > 1e-10) {
                n2 = n2/vecnorm;
            }

            n3      = arma::cross(r24,r23);
            vecnorm = arma::norm(n3);
            if (vecnorm > 1e-10) {
                n3 = n3/vecnorm;
            }

            n4      = arma::cross(r23,r13);
            vecnorm = arma::norm(n4);
            if (vecnorm > 1e-10) {
                n4 = n4/vecnorm;
            }

            Omega =  std::asin(arma::dot(n1,n2));
            Omega += std::asin(arma::dot(n2,n3));
            Omega += std::asin(arma::dot(n3,n4));
            Omega += std::asin(arma::dot(n4,n1));

            Omega = Omega*fac4pi*sgn( arma::dot( arma::cross(r34,r12),r13) );

            if (Omega==Omega) {
                (*writhe_elem)(k,l) = Omega;
                (*writhe_elem)(l,k) = Omega;
            }
        }
    }

    /*
        Calculate Gauss integral segment contributions involving the last segment
    */

    int l = num_points-1;
    j = l*point_dist;
    arma::colvec lastpos;
    if (closed_topology) {
        lastpos = bp_pos.col(0);
    }
    else {
        lastpos= bp_pos.col(num_bp-1)+triads.slice(num_bp-1).col(2)*disc_len;
    }
    for (int k=1;k<num_points-2;k++) {
        i   = k*point_dist;
        ip1 = (k+1)*point_dist;

        r12 = bp_pos.col(ip1) - bp_pos.col(i);
        r13 = bp_pos.col(j)   - bp_pos.col(i);
        r14 = lastpos         - bp_pos.col(i);
        r23 = bp_pos.col(j)   - bp_pos.col(ip1);
        r24 = lastpos         - bp_pos.col(ip1);
        r34 = lastpos         - bp_pos.col(j);

        n1      = arma::cross(r13,r14);
        vecnorm = arma::norm(n1);
        if (vecnorm > 1e-10) {
            n1 = n1/vecnorm;
        }

        n2      = arma::cross(r14,r24);
        vecnorm = arma::norm(n2);
        if (vecnorm > 1e-10) {
            n2 = n2/vecnorm;
        }

        n3      = arma::cross(r24,r23);
        vecnorm = arma::norm(n3);
        if (vecnorm > 1e-10) {
            n3 = n3/vecnorm;
        }

        n4      = arma::cross(r23,r13);
        vecnorm = arma::norm(n4);
        if (vecnorm > 1e-10) {
            n4 = n4/vecnorm;
        }

        Omega =  std::asin(arma::dot(n1,n2));
        Omega += std::asin(arma::dot(n2,n3));
        Omega += std::asin(arma::dot(n3,n4));
        Omega += std::asin(arma::dot(n4,n1));

        Omega = Omega*fac4pi*sgn( arma::dot( arma::cross(r34,r12),r13) );

        if (Omega==Omega) {
            (*writhe_elem)(k,l) = Omega;
            (*writhe_elem)(l,k) = Omega;
        }
    }

}









