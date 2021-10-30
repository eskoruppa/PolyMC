#include "Constraint.h"

Constraint::Constraint(Chain * ch,const std::vector<int> & constr_group) :
chain(ch),
constr_group(constr_group),
bp_pos(ch->get_bp_pos()),
triads(ch->get_triads()),
num_bp(ch->get_num_bp()),
num_bps(ch->get_num_bps()),
disc_len(ch->get_disc_len()),
counter(0),
counter_reject(0),
bp_pos_backup(*bp_pos),
triads_backup(*triads),
ptr_bp_pos_backup(&bp_pos_backup),
ptr_triads_backup(&triads_backup),
intrinsic_backup(true),
set_backup(true)
{
    check_valid_constr_group();
}

Constraint::Constraint(Chain * ch,const std::vector<int> & constr_group, arma::mat * backup_pos, arma::cube * backup_triads) :
chain(ch),
constr_group(constr_group),
bp_pos(ch->get_bp_pos()),
triads(ch->get_triads()),
num_bp(ch->get_num_bp()),
num_bps(ch->get_num_bps()),
disc_len(ch->get_disc_len()),
counter(0),
counter_reject(0),
ptr_bp_pos_backup(backup_pos),
ptr_triads_backup(backup_triads),
intrinsic_backup(false),
set_backup(false)
{
    check_valid_constr_group();
}

Constraint::~Constraint() {}


bool Constraint::check(const std::vector<arma::ivec>* moved) {
    counter++;
    bool accepted = check_constraint(moved);
    if (accepted) {
        if (set_backup) {
            set_current_as_backup();
        }
    }
    else {
        counter_reject++;
        revert_to_backup();
    }

//    if (accepted) {
//        std::cout << "constraint fullfilled (" << name << ")" << std::endl;
//    }
//    else {
//        std::cout << "constraint violated (" << name << ")" << std::endl;
//    }

    return accepted;
}

bool Constraint::check_constraint(const std::vector<arma::ivec>* moved) {
    /*
        This is a placeholder for the particular constraint check
    */
    return true;
}

bool Constraint::check_valid_constr_group() {
    for (unsigned i=0;i<constr_group.size();i++) {
        if (constr_group[i] < 0 || constr_group[i] >= num_bp) {
            std::cout << "Error: Invalid monomer (" <<constr_group[i] << ") selected in constraint group!" << std::endl;
            std::cout << "       The current configuration contains " << num_bp << " monomers."            << std::endl;
            return false;
        }
    }
    return true;
}

arma::mat*  Constraint::get_bp_pos_backup() {
    return &bp_pos_backup;
}

arma::cube* Constraint::get_triads_backup() {
    return &triads_backup;
}

void Constraint::set_backup_conf(arma::mat * backup_pos, arma::cube * backup_triads) {
    bp_pos_backup = *backup_pos;
    triads_backup = *backup_triads;
}

void Constraint::transfer_config(const std::vector<arma::ivec>* moved, arma::mat* bp_pos_from, arma::mat* bp_pos_to, arma::cube* triads_from, arma::cube* triads_to ) {
/*
    Copying the full arma::mat and arma::cube is much faster as counted per element as copying partial mat/cube.
    Hence if the amount of moved elements exceeds 0.2*num_bp the full config will be copied
*/
    *bp_pos_to = *bp_pos_from;
    *triads_to = *triads_from;
}

void Constraint::revert_to_backup() {
/*
    Reverts the current chain positions and triads back to the stored backup.
*/
    *bp_pos = *ptr_bp_pos_backup;
    *triads = *ptr_triads_backup;
}

void Constraint::set_current_as_backup() {
/*
    Sets the current chain positions and triads as backup.
*/
    *ptr_bp_pos_backup = *bp_pos;
    *ptr_triads_backup = *triads;
}

double Constraint::rejection_rate() {
    if (counter>0) {
        return (double)counter_reject/(double)counter;
    }
    else {
        return 0;
    }
}

void Constraint::allow_setting_backup(bool allowed) {
    set_backup = allowed;
}

const std::vector<int> Constraint::filter_double_entries(const std::vector<int> & group) {

    std::vector<int> filtered;
    for (unsigned i=0;i<group.size();i++) {
        bool already_contained = false;
        for (unsigned j=0;j<filtered.size();j++) {
            if (group[i] == filtered[j]) {
                already_contained = true;
                break;
            }
        }
        if (!already_contained) {
            filtered.push_back(group[i]);
        }
    }
    return filtered;
}

