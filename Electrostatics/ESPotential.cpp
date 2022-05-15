#include "ESPotential.h"

ESPotential::ESPotential(Chain * ch, std::vector<double> & params, double integral_dx, double rho_max) :
chain(ch),
disc_len(ch->get_disc_len()),
kT(ch->get_kT()),
params(params),
integral_dx(integral_dx),
rho_max(rho_max),
tabulate(false),
use_costheta(true),
table_N(0),
min_dist(0),
max_dist(0)
{
//    prefactor = set_prefactor();

    num_discretizations = int(std::ceil(disc_len/integral_dx));
    integral_dx = disc_len/num_discretizations;

    charge_density = CHARGE_PER_BP/DISTANCE_BETWEEN_BP;
}

ESPotential::~ESPotential() {

}

double  ESPotential::segmentpair_energy( const arma::colvec & r,
                                         const arma::colvec & t1,
                                         const arma::colvec & t2) {
    return segmentpair_energy(r,t1,t2,arma::norm(r));
}

double  ESPotential::segmentpair_energy( const arma::colvec & r,
                                         const arma::colvec & t1,
                                         const arma::colvec & t2,
                                         const double dist) {
    if (dist >= rho_max) {
        return 0;
    }
    if (dist < min_dist) {
        return 1e10;
    }

//    arma::colvec cr,ct1,ct2,ct3;
//    double cdist;
//    cr   = {4,0,0};
//    ct1  = {0,0,1};
//    ct2  = {0,1,0};
//    cdist = 4;
//
//
////    eval_double_integral( cr, ct1,ct2, cdist);
//    eval_table(cr,ct1,ct2,cdist);



    if (tabulate) {
        if (use_costheta) {
            return eval_table(r,t1,t2,dist);
        }
        else {
            return eval_table_theta(r,t1,t2,dist);
        }
    }
    return eval_double_integral(r,t1,t2,dist);
}

double ESPotential::eval_double_integral( const arma::colvec & r,
                                          const arma::colvec & t1,
                                          const arma::colvec & t2,
                                          const double dist) {

    /*
        Numerical integration of double integral over 2 segments
    */

    #ifdef DEBUG_ESPOTENTIAL
    if (std::abs(arma::norm(t1) - 1) >  1e-12) {
        throw std::invalid_argument("ESPotential::eval_double_integral(): Vector t1 not normalized!");
    }
    if (std::abs(arma::norm(t2) - 1) >  1e-12) {
        throw std::invalid_argument("ESPotential::eval_double_integral(): Vector t1 not normalized!");
    }
    #endif

    arma::colvec p,q,p1,q1,dp,dq;
    double r_norm;

    p1 =     0.5*(-disc_len + integral_dx)*t1;
    q1 = r + 0.5*(-disc_len + integral_dx)*t2;

//    p  = p1;
//    q  = q1;
    dp = t1*integral_dx;
    dq = t2*integral_dx;

    double sum = 0;

//    // 0 - 0
//    sum += integrant(arma::norm(p-q));
//
//    // 0 - j (j>0)
//    for (int j=1;j<num_discretizations;j++) {
//        q    = q+dq;
//        sum += integrant(arma::norm(p-q));
//    }
//
//    for (int i=1;i<num_discretizations;i++) {
//        p = p+dp;
//        q = q1;
//        // i (j>0) - 0
//        sum += integrant(arma::norm(p-q));
//
//        // i (j>0) - j (j>0)
//        for (int j=1;j<num_discretizations;j++) {
//            q = q+dq;
////            r    = arma::norm(p-q);
//            sum += integrant(arma::norm(p-q));
//        }
//    }
//
//    #ifdef DEBUG_ESPOTENTIAL
//    double checksum = 0;
//    for (int i=0;i<num_discretizations;i++) {
//        for (int j=0;j<num_discretizations;j++) {
//            p = p1+i*dp;
//            q = q1+j*dq;
//            checksum += integrant(arma::norm(p-q));
//        }
//    }
//    if (!equal_double(checksum*integral_dx*integral_dx,sum*integral_dx*integral_dx,1e-6)) {
//        std::cout << "sums not equal!" << std::endl;
//        std::cout << "sum:      " << sum << std::endl;
//        std::cout << "checksum: " << checksum << std::endl;
//        std::cout << "diff:     " << sum-checksum << std::endl;
//        std::exit(0);
//    }
//    #endif

    for (int i=0;i<num_discretizations;i++) {
        for (int j=0;j<num_discretizations;j++) {
            p = p1+i*dp;
            q = q1+j*dq;
            r_norm = arma::norm(p-q);
            if (r_norm<rho_max) {
                sum += integrant(r_norm);
            }
        }
    }

    sum = sum*integral_dx*integral_dx*prefactor;
    return sum;
}

double ESPotential::eval_double_integral( const arma::colvec & r,
                                          const arma::colvec & t1,
                                          const arma::colvec & t2,
                                          const double dist,
                                          double dx ) {

    /*
        Numerical integration of double integral over 2 segments
        this version allows for setting dx
    */

    int num_disc = int(std::ceil(disc_len/dx));
    dx       = disc_len/num_disc;

    #ifdef DEBUG_ESPOTENTIAL
    if (std::abs(arma::norm(t1) - 1) >  1e-12) {
        throw std::invalid_argument("ESPotential::eval_double_integral(): Vector t1 not normalized!");
    }
    if (std::abs(arma::norm(t2) - 1) >  1e-12) {
        throw std::invalid_argument("ESPotential::eval_double_integral(): Vector t1 not normalized!");
    }
    #endif

    arma::colvec p,q,p1,q1,dp,dq;

    p1 =     0.5*(-disc_len + dx)*t1;
    q1 = r + 0.5*(-disc_len + dx)*t2;

    dp = t1*dx;
    dq = t2*dx;

    double sum = 0;

    for (int i=0;i<num_disc;i++) {
        for (int j=0;j<num_disc;j++) {
            p = p1+i*dp;
            q = q1+j*dq;
            sum += integrant(arma::norm(p-q));
        }
    }

    sum = sum*dx*dx*prefactor;
    return sum;
}


std::string ESPotential::get_type() {
    return type;
}

double ESPotential::set_prefactor() {
    throw "ESPotential::set_prefactor(): The base class ESPotential::set_prefactor should never be called.";
}
double ESPotential::integrant(double r) {
    /*
        Define the integrant
    */
    throw "ESPotential::integrant(): The base class ESPotential::integrant should never be called.";
}


bool    ESPotential::tabulate_interactions(int table_N,double min_dist,double max_dist,double memory_limit,bool use_costheta) {
    if (table_N==0) {
        return false;
    }
    tabulate     = true;
    this->use_costheta  = use_costheta;
    this->table_N       = table_N;
    this->min_dist      = min_dist;
    this->max_dist      = max_dist;
    this->memory_limit  = memory_limit;
    if (use_costheta) {
        init_tabulation_costheta();
    }
    else {
        init_tabulation_theta();
    }
    return true;
}


///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////// TABULATE DATA /////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void ESPotential::init_tabulation_costheta() {

    std::cout << std::endl;
    std::cout << "######################################"   << std::endl;
    std::cout << "#### INIT ELECTROSTATICS TABLE #######"   << std::endl;
    std::cout << "######################################"   << std::endl;
    int max_N = int(std::pow(1e6*memory_limit*1./8,0.25));
    if (table_N > max_N) {
        std::cout << " Warning: Current number of table entries would exceed maximum memory usage." << std::endl;
        std::cout << "          -> setting number of entries to " << max_N << std::endl;
        table_N = max_N;
    }
    double memory_use = 8.*table_N*table_N*table_N*table_N/1000000;
    std::cout << "   type:              " << type   << std::endl;
    std::cout << "   integral_dx:       " << integral_dx   << std::endl;
    std::cout << "   entries per dim:   " << table_N   << std::endl;
    std::cout << "   use cos theta:     " << use_costheta  << std::endl;
    std::cout << "   min distance:      " << min_dist  << std::endl;
    std::cout << "   max distance:      " << max_dist  << std::endl;
    std::cout << "   approx memory use: " << memory_use << " MB"   << std::endl;
    std::cout << "######################################"   << std::endl;
    std::cout << std::endl;

    table_intervals     = table_N-1;

    table_dist_range    = max_dist-min_dist;
    table_dist_step     = table_dist_range/(table_intervals);
    table_dist_vals     = arma::zeros(table_N);

    table_costheta_step = (2.)/(table_intervals);
    table_costheta_vals = arma::zeros(table_N);

    for (int i=0;i<table_N-1;i++) {
        table_dist_vals[i]     = min_dist+i*table_dist_step;
        table_costheta_vals[i] = -1 + i*table_costheta_step;
    }
    table_dist_vals[table_N-1]     = max_dist;
    table_costheta_vals[table_N-1] = 1;

    // init table
    for (int i=0;i<table_N;i++) {
        arma::cube new_entry = arma::zeros(table_N,table_N,table_N);
        table.push_back(new_entry);
    }

    double cth1,cth2,cth3;
    double th1,th2,th3;
    double dist;

    double value;

    arma::mat R1,R2,R3;
    arma::colvec t1,t2;
    arma::colvec t1_init,t2_init;
    arma::colvec t1_xrot,t2_xrot;
    arma::colvec t1_twist;

    t1_init = {1,0,0};
    t2_init = {1,0,0};

    arma::colvec r;
    arma::colvec rho;
    r   = {1,0,0};
    rho = {1,0,0};

    // Progress Printing
    long long next    = 2;
    long long counter = 0;
    long long total   = table_N*table_N*table_N*table_N;

    std::cout << " Building lookup table for " << total << " segment integrals ..." << std::endl;
    std::cout << "           |0%____________________________________________100%|" << std::endl;
    std::cout << " Progress:  " << std::flush;
    /////////////////////////////////////////


    for (int i=0;i<table_N;i++) {
        cth1 = -1+i*table_costheta_step;
        th1  = std::acos( cth1 );
        R1   = Rotz(th1);
        t1_xrot = R1*t1_init;


        for (int j=0;j<table_N;j++) {
            cth2 = -1+j*table_costheta_step;
            th2  = std::acos( cth2);
            R2   = Rotz(th2);
            t2   = R2*t2_init;

            for (int k=0;k<table_N;k++) {
                cth3 = -1+k*table_costheta_step;
                th3  = std::acos( cth3 );
                R3   = Rotx(th3);

                t1 = R3*t1_xrot;

//                std::cout << "############" << std::endl;

                #ifdef DEBUG_ESPOTENTIAL
                ////////////////////////
                // checks
                double check_cth1 = arma::dot(rho,t1);
                double check_cth2 = arma::dot(rho,t2);
                double check_cth3;

                arma::colvec cross1 = arma::cross(t1,rho);
                arma::colvec cross2 = arma::cross(t2,rho);
                double nc1 = arma::norm(cross1);
                double nc2 = arma::norm(cross2);
                if (nc1<1e-14 || nc2<1e-14) {
                    check_cth3 = 0;
                }
                else {
                    check_cth3 = arma::dot(cross1,cross2)/nc1/nc2;
                }

                if (std::abs(check_cth1-cth1) > 1e-12) {
                    std::cout << "cth1 = " << cth1 << " - " << check_cth1 <<  std::endl;
                    std::cout << "wrongly assigned cth1" << std::endl;
                    std::exit(0);
                }
                if (std::abs(check_cth1-cth1) > 1e-12) {
                    std::cout << "cth2 = " << cth2 << " - " << check_cth2 <<  std::endl;
                    std::cout << "wrongly assigned cth2" << std::endl;
                    std::exit(0);
                }
                if (std::abs(check_cth3-cth3) > 1e-12 && std::abs(std::abs(cth1)-1) > 1e-12 && std::abs(std::abs(cth2)-1) > 1e-12) {
                    std::cout << "----------" << std::endl;
                    std::cout << "cth3 = " << cth3 << " - " << check_cth3 <<  std::endl;
                    std::cout << "cth1 = " << cth1 <<  std::endl;
                    std::cout << "cth2 = " << cth2 <<  std::endl;
                    std::cout << "wrongly assigned cth3" << std::endl;
                    std::exit(0);
                }
                #endif

                #ifdef DEBUG_ESPOTENTIAL
                if (!equal_double(arma::norm(t1),1)) {
                    std::cout << "t1 not normed" << std::endl;
                    std::exit(0);
                }
                if (!equal_double(arma::norm(t2),1)) {
                    std::cout << "t2 not normed" << std::endl;
                    std::exit(0);
                }
                #endif

                for (int d=0;d<table_N;d++) {
                    dist  = min_dist+d*table_dist_step;
                    r[0] = dist;

                    ////////////////////////////prefactor*
                    value = eval_double_integral( r, t1,t2, dist);
                    set_table_val(value,d,i,j,k);

                    /////////////////////////////////////////
                    counter++;
                    if (counter>= 1.0*next/100.*total) {
                        std::cout << "\u2588" << std::flush;
                        next += 2;
                    }
                    /////////////////////////////////////////
                }
            }
        }
    }
    /////////////////////////////////////////
    std::cout << std::endl << "... done!" << std::endl;
    std::cout << std::endl;

}


double ESPotential::eval_table( const arma::colvec & r,
                                const arma::colvec & t1,
                                const arma::colvec & t2,
                                const double dist) {

    double cth1,cth2,cth3;
    arma::colvec rho;
    arma::colvec cross1,cross2;
    double nc1,nc2;
    int id1,id2,id3,idd;
    int id1p1,id2p1,id3p1,iddp1;

    rho = r/arma::norm(r);
    cth1 = arma::dot(rho,t1);
    cth2 = arma::dot(rho,t2);

    cross1 = arma::cross(t1,rho);
    cross2 = arma::cross(t2,rho);
    nc1 = arma::norm(cross1);
    nc2 = arma::norm(cross2);
    if (nc1<ZERO_CROSS || nc2<ZERO_CROSS) {
        cth3 = 0;
    }
    else {
        cth3 = arma::dot(cross1,cross2)/nc1/nc2;
    }

    idd = (dist-min_dist)/table_dist_step;
    id1 = (cth1+1)/table_costheta_step;
    id2 = (cth2+1)/table_costheta_step;
    id3 = (cth3+1)/table_costheta_step;


//    if (id1<0) {
//        id1  = 0;
//        cth1 = table_costheta_vals(0);
//    }
//    if (id2<0) {
//        id2  = 0;
//        cth2 = table_costheta_vals(0);
//    }
//    if (id3<0) {
//        id3  = 0;
//        cth3 = table_costheta_vals(0);
//    }


    iddp1   = idd+1;
    id1p1   = id1+1;
    id2p1   = id2+1;
    id3p1   = id3+1;

    #ifdef DEBUG_ESPOTENTIAL
    if (table_dist_vals(idd) > dist || table_dist_vals(iddp1) < dist) {
        std::cout <<  "wrongly assigned distance id" << std::endl;
        std::exit(0);
    }
    if (table_costheta_vals(id1) > cth1 || table_dist_vals(id1p1) < cth1) {
        if (!equal_double(-1,cth1,1e-10)) {
            std::cout <<  "wrongly assigned cth1 id" << std::endl;
            std::cout <<  cth1 << std::endl;
            std::cout <<  table_costheta_vals(id1) << std::endl;
            std::cout <<  table_costheta_vals(id1p1) << std::endl;
            std::exit(0);
        }
    }
    if (table_costheta_vals(id2) > cth2 || table_dist_vals(id2p1) < cth2) {
        if (!equal_double(-1,cth2,1e-10)) {
            std::cout <<  "wrongly assigned cth2 id" << std::endl;
            std::cout <<  cth2 << std::endl;
            std::cout <<  table_costheta_vals(id2) << std::endl;
            std::cout <<  table_costheta_vals(id2p1) << std::endl;
            std::exit(0);
        }
    }
    if (table_costheta_vals(id3) > cth3 || table_dist_vals(id3p1) < cth3) {
        if (!equal_double(-1,cth3,1e-10)) {
            std::cout <<  "wrongly assigned cth3 id" << std::endl;
            std::cout <<  cth3 << std::endl;
            std::cout <<  table_costheta_vals(id3) << std::endl;
            std::cout <<  table_costheta_vals(id3p1) << std::endl;

            std::cout <<  table_costheta_vals(0)-cth3 << std::endl;
            std::exit(0);
        }
    }
    #endif

    double ld,lcth1,lcth2,lcth3;
    double dd_lo,ct1_lo,ct2_lo,ct3_lo;
    double dd_up,ct1_up,ct2_up,ct3_up;

    ld      = table_dist_vals(idd);
    lcth1   = table_costheta_vals(id1);
    lcth2   = table_costheta_vals(id2);
    lcth3   = table_costheta_vals(id3);


    dd_lo  = (dist-ld)/table_dist_step;
    ct1_lo = (cth1-lcth1)/table_costheta_step;
    ct2_lo = (cth2-lcth2)/table_costheta_step;
    ct3_lo = (cth3-lcth3)/table_costheta_step;


    #ifdef DEBUG_ESPOTENTIAL
    if (dd_lo < 0 || dd_lo > 1) {
        std::cout <<  "dist fraction out of bounds" << std::endl;
        std::exit(0);
    }
    if (ct1_lo < 0 || ct1_lo > 1) {
        if (!equal_double(-1,cth1,1e-10)) {
            std::cout <<  "costheta1 fraction out of bounds" << std::endl;
            std::cout <<  ct1_lo << std::endl;
            std::exit(0);
        }
    }
    if (ct2_lo < 0 || ct2_lo > 1) {
        if (!equal_double(-1,cth2,1e-10)) {
            std::cout <<  "costheta2 fraction out of bounds" << std::endl;
            std::cout <<  ct2_lo << std::endl;
            std::exit(0);
        }
    }
    if (ct3_lo < 0 || ct3_lo > 1) {
        if (!equal_double(-1,cth3,1e-10)) {
            std::cout <<  "costheta3 fraction out of bounds" << std::endl;
            std::cout <<  ct3_lo << std::endl;

            std::cout <<  cth3 << std::endl;
            std::cout <<  table_costheta_vals(id3) << std::endl;
            std::cout <<  table_costheta_vals(id3p1) << std::endl;
            std::cout <<  id3 << std::endl;

            std::cout <<  table_costheta_vals(0)-cth3 << std::endl;
            std::exit(0);
        }
    }
    #endif

    dd_up  = 1-dd_lo;
    ct1_up = 1-ct1_lo;
    ct2_up = 1-ct2_lo;
    ct3_up = 1-ct3_lo;

    double value = 0;
    value += dd_lo * ct1_lo * ct2_lo * ct3_lo * get_table_val(idd, id1, id2, id3)  ;
    value += dd_lo * ct1_lo * ct2_lo * ct3_up * get_table_val(idd, id1, id2, id3p1);

    value += dd_lo * ct1_lo * ct2_up * ct3_lo * get_table_val(idd, id1, id2p1, id3)  ;
    value += dd_lo * ct1_lo * ct2_up * ct3_up * get_table_val(idd, id1, id2p1, id3p1);

    value += dd_lo * ct1_up * ct2_lo * ct3_lo * get_table_val(idd, id1p1, id2, id3)  ;
    value += dd_lo * ct1_up * ct2_lo * ct3_up * get_table_val(idd, id1p1, id2, id3p1);

    value += dd_lo * ct1_up * ct2_up * ct3_lo * get_table_val(idd, id1p1, id2p1, id3)  ;
    value += dd_lo * ct1_up * ct2_up * ct3_up * get_table_val(idd, id1p1, id2p1, id3p1);


    value += dd_up * ct1_lo * ct2_lo * ct3_lo * get_table_val(iddp1, id1, id2, id3)  ;
    value += dd_up * ct1_lo * ct2_lo * ct3_up * get_table_val(iddp1, id1, id2, id3p1);

    value += dd_up * ct1_lo * ct2_up * ct3_lo * get_table_val(iddp1, id1, id2p1, id3)  ;
    value += dd_up * ct1_lo * ct2_up * ct3_up * get_table_val(iddp1, id1, id2p1, id3p1);

    value += dd_up * ct1_up * ct2_lo * ct3_lo * get_table_val(iddp1, id1p1, id2, id3)  ;
    value += dd_up * ct1_up * ct2_lo * ct3_up * get_table_val(iddp1, id1p1, id2, id3p1);

    value += dd_up * ct1_up * ct2_up * ct3_lo * get_table_val(iddp1, id1p1, id2p1, id3)  ;
    value += dd_up * ct1_up * ct2_up * ct3_up * get_table_val(iddp1, id1p1, id2p1, id3p1);


//    double check_value = eval_double_integral( r, t1,t2, dist);
//    double check_few = eval_double_integral( r, t1,t2, dist,1);
//    std::cout << "\n\n" << idd << " " << id1 << " " << id2 << " " << id3 << std::endl;
//    std::cout << "value     = " << value << std::endl;
//    std::cout << "sval      = " << get_table_val(idd, id1, id2, id3) << std::endl;
//    std::cout << "check     = " << check_value << std::endl;
//    std::cout << "check_few = " << check_few << std::endl;

//    for (int i=0;i<2;i++) {
//        for (int j=0;j<2;j++) {
//            for (int k=0;k<2;k++) {
//                for (int l=0;l<2;l++) {
//                    std::cout << " vals " << get_table_val(idd+i, id1+j, id2+k, id3+l) << std::endl;
//                }
//            }
//        }
//    }

    return value;
}








//////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////


void ESPotential::init_tabulation_theta() {

    std::cout << std::endl;
    std::cout << "######################################"   << std::endl;
    std::cout << "#### INIT ELECTROSTATICS TABLE #######"   << std::endl;
    std::cout << "######################################"   << std::endl;
    int max_N = int(std::pow(1e6*memory_limit*1./8,0.25));
    if (table_N > max_N) {
        std::cout << " Warning: Current number of table entries would exceed maximum memory usage." << std::endl;
        std::cout << "          -> setting number of entries to " << max_N << std::endl;
        table_N = max_N;
    }
    double memory_use = 8.*table_N*table_N*table_N*table_N/1000000;
    std::cout << "   type:              " << type   << std::endl;
    std::cout << "   integral_dx:       " << integral_dx   << std::endl;
    std::cout << "   entries per dim:   " << table_N   << std::endl;
    std::cout << "   use cos theta:     " << use_costheta  << std::endl;
    std::cout << "   min distance:      " << min_dist  << std::endl;
    std::cout << "   max distance:      " << max_dist  << std::endl;
    std::cout << "   approx memory use: " << memory_use << " MB"   << std::endl;
    std::cout << "######################################"   << std::endl;
    std::cout << std::endl;

    table_intervals     = table_N-1;

    table_dist_range    = max_dist-min_dist;
    table_dist_step     = table_dist_range/(table_intervals);
    table_dist_vals     = arma::zeros(table_N);

    table_costheta_step = (M_PI)/(table_intervals);
    table_costheta_vals = arma::zeros(table_N);

    for (int i=0;i<table_N-1;i++) {
        table_dist_vals[i]     = min_dist+i*table_dist_step;
        table_costheta_vals[i] = i*table_costheta_step;
    }
    table_dist_vals[table_N-1]     = max_dist;
    table_costheta_vals[table_N-1] = M_PI;

    // init table
    for (int i=0;i<table_N;i++) {
        arma::cube new_entry = arma::zeros(table_N,table_N,table_N);
        table.push_back(new_entry);
    }

    double cth1,cth2,cth3;
    double th1,th2,th3;
    double dist;

    double value;

    arma::mat R1,R2,R3;
    arma::colvec t1,t2;
    arma::colvec t1_init,t2_init;
    arma::colvec t1_xrot,t2_xrot;
    arma::colvec t1_twist;

    t1_init = {1,0,0};
    t2_init = {1,0,0};

    arma::colvec r;
    arma::colvec rho;
    r   = {1,0,0};
    rho = {1,0,0};

    // Progress Printing
    long long next    = 2;
    long long counter = 0;
    long long total   = table_N*table_N*table_N*table_N;

    std::cout << " Building lookup table for " << total << " segment integrals ..." << std::endl;
    std::cout << "           |0%____________________________________________100%|" << std::endl;
    std::cout << " Progress:  " << std::flush;
    /////////////////////////////////////////


    for (int i=0;i<table_N;i++) {
        th1  = i*table_costheta_step;
        cth1 = std::cos(th1);
        R1   = Rotz(th1);
        t1_xrot = R1*t1_init;


        for (int j=0;j<table_N;j++) {
            th2  = j*table_costheta_step;
            cth2 = std::cos(th2);
            R2   = Rotz(th2);
            t2   = R2*t2_init;

            for (int k=0;k<table_N;k++) {
                th3  = k*table_costheta_step;
                cth3 = std::cos(th3);
                R3   = Rotx(th3);

                t1 = R3*t1_xrot;

//                std::cout << "############" << std::endl;

                #ifdef DEBUG_ESPOTENTIAL
                ////////////////////////
                // checks
                double check_cth1 = arma::dot(rho,t1);
                double check_cth2 = arma::dot(rho,t2);
                double check_cth3;
                double check_th1,check_th2,check_th3;

                arma::colvec cross1 = arma::cross(t1,rho);
                arma::colvec cross2 = arma::cross(t2,rho);
                double nc1 = arma::norm(cross1);
                double nc2 = arma::norm(cross2);
                if (nc1<1e-14 || nc2<1e-14) {
                    check_cth3 = 0;
                }
                else {
                    check_cth3 = arma::dot(cross1,cross2)/nc1/nc2;
                }
                check_th1 = std::acos(check_cth1);
                check_th2 = std::acos(check_cth2);
                check_th3 = std::acos(check_cth3);

                if (std::abs(check_th1-th1) > 1e-12) {
                    std::cout << "th1 = " << th1 << " - " << check_th1 <<  std::endl;
                    std::cout << "wrongly assigned th1" << std::endl;
                    std::exit(0);
                }
                if (std::abs(check_th1-th1) > 1e-12) {
                    std::cout << "th2 = " << th2 << " - " << check_th2 <<  std::endl;
                    std::cout << "wrongly assigned th2" << std::endl;
                    std::exit(0);
                }
                if (std::abs(check_th3-th3) > 1e-7 && std::abs(std::abs(cth1)-1) > 1e-12 && std::abs(std::abs(cth2)-1) > 1e-12) {
                    std::cout << "----------" << std::endl;
                    std::cout << "th3 = " << th3 << " - " << check_th3 <<  std::endl;
                    std::cout << "cth1 = " << cth1 <<  std::endl;
                    std::cout << "cth2 = " << cth2 <<  std::endl;
                    std::cout << "wrongly assigned th3" << std::endl;
                    std::exit(0);
                }
                #endif

                #ifdef DEBUG_ESPOTENTIAL
                if (!equal_double(arma::norm(t1),1)) {
                    std::cout << "t1 not normed" << std::endl;
                    std::exit(0);
                }
                if (!equal_double(arma::norm(t2),1)) {
                    std::cout << "t2 not normed" << std::endl;
                    std::exit(0);
                }
                #endif

                for (int d=0;d<table_N;d++) {
                    dist  = min_dist+d*table_dist_step;
                    r[0] = dist;

                    ////////////////////////////prefactor*
                    value = eval_double_integral( r, t1,t2, dist);
                    set_table_val(value,d,i,j,k);

                    /////////////////////////////////////////
                    counter++;
                    if (counter>= 1.0*next/100.*total) {
                        std::cout << "\u2588" << std::flush;
                        next += 2;
                    }
                    /////////////////////////////////////////
                }
            }
        }
    }
    /////////////////////////////////////////
    std::cout << std::endl << "... done!" << std::endl;
    std::cout << std::endl;
}


double ESPotential::eval_table_theta(   const arma::colvec & r,
                                        const arma::colvec & t1,
                                        const arma::colvec & t2,
                                        const double dist) {

    double th1,th2,th3;
    double cth1,cth2,cth3;
    arma::colvec rho;
    arma::colvec cross1,cross2;
    double nc1,nc2;
    int id1,id2,id3,idd;
    int id1p1,id2p1,id3p1,iddp1;

    rho = r/arma::norm(r);
    cth1 = arma::dot(rho,t1);
    cth2 = arma::dot(rho,t2);

    cross1 = arma::cross(t1,rho);
    cross2 = arma::cross(t2,rho);
    nc1 = arma::norm(cross1);
    nc2 = arma::norm(cross2);
    if (nc1<ZERO_CROSS || nc2<ZERO_CROSS) {
        cth3 = 0;
    }
    else {
        cth3 = arma::dot(cross1,cross2)/nc1/nc2;
    }
    if (cth1<-1) {
        cth1 = -1+1e-6;
    }
    if (cth1>1) {
        cth1 = 1-1e-6;
    }
    if (cth2<-1) {
        cth2 = -1+1e-6;
    }
    if (cth2>1) {
        cth2 = 1-1e-6;
    }
    if (cth3<-1) {
        cth3 = -1+1e-6;
    }
    if (cth3>1) {
        cth3 = 1-1e-6;
    }

    th1 = std::acos(cth1);
    th2 = std::acos(cth2);
    th3 = std::acos(cth3);


    idd = (dist-min_dist)/table_dist_step;
    id1 = (th1)/table_costheta_step;
    id2 = (th2)/table_costheta_step;
    id3 = (th3)/table_costheta_step;

    int num_intervals = table_N-1;
    if (id1 >= num_intervals) {
        id1 = num_intervals-1;
    }
    if (id2 >= num_intervals) {
        id2 = num_intervals-1;
    }
    if (id3 >= num_intervals) {
        id3 = num_intervals-1;
    }


//    if (id1<0) {
//        id1  = 0;
//        cth1 = table_costheta_vals(0);
//    }
//    if (id2<0) {
//        id2  = 0;
//        cth2 = table_costheta_vals(0);
//    }
//    if (id3<0) {
//        id3  = 0;
//        cth3 = table_costheta_vals(0);
//    }


    iddp1   = idd+1;
    id1p1   = id1+1;
    id2p1   = id2+1;
    id3p1   = id3+1;

//    std::cout <<  "id" << std::endl;
//    std::cout <<  id1 << std::endl;
//    std::cout <<  id2 << std::endl;
//    std::cout <<  id3 << std::endl;
//
//    std::cout <<  "th" << std::endl;
//    std::cout <<  th1 << std::endl;
//    std::cout <<  th2 << std::endl;
//    std::cout <<  th3 << std::endl;
//
//    std::cout <<  "cth" << std::endl;
//    std::cout <<  cth1 << std::endl;
//    std::cout <<  cth2 << std::endl;
//    std::cout <<  cth3 << std::endl;

    #ifdef DEBUG_ESPOTENTIAL
    if (table_dist_vals(idd) > dist || table_dist_vals(iddp1) < dist) {
        std::cout <<  "wrongly assigned distance id" << std::endl;
        std::exit(0);
    }
    if (table_costheta_vals(id1) > th1 || table_dist_vals(id1p1) < th1) {
//        if (!equal_double(-1,cth1,1e-10)) {
            std::cout <<  "wrongly assigned th1 id" << std::endl;
            std::cout <<  th1 << std::endl;
            std::cout <<  table_costheta_vals(id1) << std::endl;
            std::cout <<  table_costheta_vals(id1p1) << std::endl;
            std::exit(0);
//        }
    }
    if (table_costheta_vals(id2) > th2 || table_dist_vals(id2p1) < th2) {
//        if (!equal_double(-1,cth2,1e-10)) {
            std::cout <<  "wrongly assigned th2 id" << std::endl;
            std::cout <<  th2 << std::endl;
            std::cout <<  table_costheta_vals(id2) << std::endl;
            std::cout <<  table_costheta_vals(id2p1) << std::endl;
            std::exit(0);
//        }
    }
    if (table_costheta_vals(id3) > th3 || table_dist_vals(id3p1) < th3) {
//        if (!equal_double(-1,cth3,1e-10)) {
            std::cout <<  "wrongly assigned th3 id" << std::endl;
            std::cout <<  th3 << std::endl;
            std::cout <<  table_costheta_vals(id3) << std::endl;
            std::cout <<  table_costheta_vals(id3p1) << std::endl;

            std::cout <<  table_costheta_vals(0)-th3 << std::endl;
            std::exit(0);
//        }
    }
    #endif

    double ld,lcth1,lcth2,lcth3;
    double dd_lo,ct1_lo,ct2_lo,ct3_lo;
    double dd_up,ct1_up,ct2_up,ct3_up;

    ld      = table_dist_vals(idd);
    lcth1   = table_costheta_vals(id1);
    lcth2   = table_costheta_vals(id2);
    lcth3   = table_costheta_vals(id3);


    dd_lo  = (dist-ld)/table_dist_step;
    ct1_lo = (th1-lcth1)/table_costheta_step;
    ct2_lo = (th2-lcth2)/table_costheta_step;
    ct3_lo = (th3-lcth3)/table_costheta_step;


    #ifdef DEBUG_ESPOTENTIAL
    if (dd_lo < 0 || dd_lo > 1) {
        std::cout <<  "dist fraction out of bounds" << std::endl;
        std::exit(0);
    }
    if (ct1_lo < 0 || ct1_lo > 1) {
//        if (!equal_double(-1,cth1,1e-10)) {
            std::cout <<  "costheta1 fraction out of bounds" << std::endl;
            std::cout <<  ct1_lo << std::endl;
            std::exit(0);
//        }
    }
    if (ct2_lo < 0 || ct2_lo > 1) {
//        if (!equal_double(-1,cth2,1e-10)) {
            std::cout <<  "costheta2 fraction out of bounds" << std::endl;
            std::cout <<  ct2_lo << std::endl;
            std::exit(0);
//        }
    }
    if (ct3_lo < 0 || ct3_lo > 1) {
//        if (!equal_double(-1,cth3,1e-10)) {
            std::cout <<  "costheta3 fraction out of bounds" << std::endl;
            std::cout <<  ct3_lo << std::endl;

            std::cout <<  cth3 << std::endl;
            std::cout <<  table_costheta_vals(id3) << std::endl;
            std::cout <<  table_costheta_vals(id3p1) << std::endl;
            std::cout <<  id3 << std::endl;

            std::cout <<  table_costheta_vals(0)-cth3 << std::endl;
            std::exit(0);
//        }
    }
    #endif

    dd_up  = 1-dd_lo;
    ct1_up = 1-ct1_lo;
    ct2_up = 1-ct2_lo;
    ct3_up = 1-ct3_lo;

    double value = 0;

    value += dd_lo * ct1_lo * ct2_lo * ct3_lo * get_table_val(idd, id1, id2, id3)  ;
    value += dd_lo * ct1_lo * ct2_lo * ct3_up * get_table_val(idd, id1, id2, id3p1);

    value += dd_lo * ct1_lo * ct2_up * ct3_lo * get_table_val(idd, id1, id2p1, id3)  ;
    value += dd_lo * ct1_lo * ct2_up * ct3_up * get_table_val(idd, id1, id2p1, id3p1);

    value += dd_lo * ct1_up * ct2_lo * ct3_lo * get_table_val(idd, id1p1, id2, id3)  ;
    value += dd_lo * ct1_up * ct2_lo * ct3_up * get_table_val(idd, id1p1, id2, id3p1);

    value += dd_lo * ct1_up * ct2_up * ct3_lo * get_table_val(idd, id1p1, id2p1, id3)  ;
    value += dd_lo * ct1_up * ct2_up * ct3_up * get_table_val(idd, id1p1, id2p1, id3p1);


    value += dd_up * ct1_lo * ct2_lo * ct3_lo * get_table_val(iddp1, id1, id2, id3)  ;
    value += dd_up * ct1_lo * ct2_lo * ct3_up * get_table_val(iddp1, id1, id2, id3p1);

    value += dd_up * ct1_lo * ct2_up * ct3_lo * get_table_val(iddp1, id1, id2p1, id3)  ;
    value += dd_up * ct1_lo * ct2_up * ct3_up * get_table_val(iddp1, id1, id2p1, id3p1);

    value += dd_up * ct1_up * ct2_lo * ct3_lo * get_table_val(iddp1, id1p1, id2, id3)  ;
    value += dd_up * ct1_up * ct2_lo * ct3_up * get_table_val(iddp1, id1p1, id2, id3p1);

    value += dd_up * ct1_up * ct2_up * ct3_lo * get_table_val(iddp1, id1p1, id2p1, id3)  ;
    value += dd_up * ct1_up * ct2_up * ct3_up * get_table_val(iddp1, id1p1, id2p1, id3p1);


//    double check_value = eval_double_integral( r, t1,t2, dist);
//    double check_few = eval_double_integral( r, t1,t2, dist,1);
//    std::cout << "\n\n" << idd << " " << id1 << " " << id2 << " " << id3 << std::endl;
//    std::cout << "value     = " << value << std::endl;
////    std::cout << "sval      = " << get_table_val(idd, id1, id2, id3) << std::endl;
//    std::cout << "check     = " << check_value << std::endl;
////    std::cout << "check_few = " << check_few << std::endl;
//
//    std::cout << "frac = " << (value-check_value)/check_value*100 << " %" << std::endl;


//    for (int i=0;i<2;i++) {
//        for (int j=0;j<2;j++) {
//            for (int k=0;k<2;k++) {
//                for (int l=0;l<2;l++) {
//                    std::cout << " vals " << get_table_val(idd+i, id1+j, id2+k, id3+l) << std::endl;
//                }
//            }
//        }
//    }

//    std::exit(0);

    return value;
}


