#include "Dump_Closure.h"

Dump_Closure::Dump_Closure(        
        Chain * ch, 
        int N_dump, 
        const std::string& filename, 
        double dist_threshold,
        double angle_threshold,
        double twist_threshold
        )
: Dump(ch,N_dump,filename,false), 
use_dist_criterion(false),
use_angle_criterion(false),
use_twist_criterion(false),
dist_threshold(dist_threshold), 
angle_threshold(angle_threshold), 
twist_threshold(twist_threshold)
{
    if (dist_threshold > 0) {
        use_dist_criterion = true;
    }
    if (angle_threshold > 0) {
        use_angle_criterion = true;
    }
    if (twist_threshold > 0) {
        use_twist_criterion = true;
    }
    if (!use_dist_criterion && !use_angle_criterion && !use_twist_criterion) {
        active = false;
    }
    else{
        active = true;
        init_store();
        if (!app) {
            std::ofstream ofstr;
            ofstr.open(fn, std::ofstream::out | std::ofstream::trunc);
            ofstr.close();
        }
    }

}
Dump_Closure::~Dump_Closure() {}

void Dump_Closure::prod_dump() {
    if (!active) {return;}

    // calculate distance
    double dist = arma::norm(chain->get_bp_pos()->col(num_bp-1)-chain->get_bp_pos()->col(0));
    double angle = std::acos(arma::dot(chain->get_triads()->slice(0).col(2),chain->get_triads()->slice(num_bp-1).col(2)));
    arma::vec Omega = ExtractTheta(chain->get_triads()->slice(num_bp-1).t() * chain->get_triads()->slice(0) );
    double twist = Omega(2);
    if ( (!use_dist_criterion || dist <= dist_threshold) 
        && (!use_angle_criterion || angle <= angle_threshold) 
        && (!use_twist_criterion || std::abs(twist) <= twist_threshold) 
        ) 
    {
        stored_dumps(store_counter,0) = dist;
        stored_dumps(store_counter,1) = angle;
        stored_dumps(store_counter,2) = twist;
        store_counter++;
        if (store_counter >= CLOSURE_NUM_STORES) {
            write2file();
        }
    }
}

void Dump_Closure::final_dump() {
    write2file();
}

void Dump_Closure::init_store() {
    store_counter = 0;
    stored_dumps = arma::zeros(CLOSURE_NUM_STORES,3);
}

void Dump_Closure::write2file() {
    std::ofstream ofstr;
    ofstr.open(fn, std::ofstream::out | std::ofstream::app);

    for (unsigned i=0;i<store_counter;i++) {
        ofstr << stored_dumps(i,0) << " " << stored_dumps(i,1) << " " << stored_dumps(i,2) << "\n";
    }
    ofstr.close();
    store_counter = 0;
}

bool Dump_Closure::is_active() {
    return active;
}

