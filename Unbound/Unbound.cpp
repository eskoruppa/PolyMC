#include "Unbound.h"


Unbound::Unbound(Chain * ch) :
chain(ch)
{
    std::vector<int> ftm(chain->get_num_bp(),default_type);
    full_type_map = ftm;
}

Unbound::~Unbound() {
    for (unsigned i=0;i<type_groups.size();i++) {
        delete type_groups[i];
    }
    for (unsigned i=0;i<pairs.size();i++) {
        delete pairs[i];
    }
}

void Unbound::set_backup(){
    backup_full_type_map = full_type_map;
}

void Unbound::revert_to_backup() {
    full_type_map = backup_full_type_map;
    assign_type_group_elements();
}

/*---------------------------------------------------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------------------------------------------------*/


void Unbound::submit_type_changes(const std::vector<std::vector<int>> & changes) {
    #ifdef UNBOUND_DEBUG
    if (type_change_acceptance_pending) {
        std::cout << "Error attempting to submit type changes while type change acceptance is still pending!" << std::endl;
        std::exit(0);
    }
    #endif
    for (unsigned i=0;i<changes.size();i++) {
        type_changes.push_back(changes[i]);
    }
}

double Unbound::eval_type_changes(const arma::mat * pos) {
    assign_trial_elements();
    type_change_acceptance_pending = true;

    double Delta_E = 0;
    for (unsigned p=0;p<pairs.size();p++) {
        Delta_E += pairs[p]->Eval_change_type(pos);
    }
    return Delta_E;
}

void Unbound::set_type_changes(bool accept) {
    type_change_acceptance_pending = false;
    if (accept) {
        /*
            Shift elements_trial to elements in type_groups
        */
        for (unsigned tg=0;tg<type_groups.size();tg++){
            type_groups[tg]->elements = type_groups[tg]->elements_trial;
        }
        for (unsigned tc=0;tc<type_changes.size();tc++) {
            full_type_map[type_changes[tc][TYPE_CHANGE_MONOMER_ID]] = type_changes[tc][TYPE_CHANGE_CHANGE_TO];
        }
    }
    type_changes.clear();
    type_changes.shrink_to_fit();
}



void Unbound::assign_trial_elements() {
#ifdef UNBOUND_DEBUG
    /*
        Check for repeated atom changes
    */
    for (unsigned i=0;i<type_changes.size()-1;i++) {
        for (unsigned j=i+1;j<type_changes.size();j++) {
            if (type_changes[i][TYPE_CHANGE_MONOMER_ID] == type_changes[j][TYPE_CHANGE_MONOMER_ID]) {
                std::cout << "Error: Unbound::eval_type_changes(): atom type change was submitted twice for monomer " << type_changes[j][TYPE_CHANGE_MONOMER_ID] << std::endl;
                std::exit(0);
            }
        }
    }
    #endif

    /*
        Sort type_changes
    */
    std::vector<int> tmp;
    for (unsigned i=0;i<type_changes.size()-1;i++) {
        for (unsigned j=i+1;j<type_changes.size();j++) {
            if (type_changes[i][TYPE_CHANGE_MONOMER_ID] > type_changes[j][TYPE_CHANGE_MONOMER_ID]) {
                tmp = type_changes[i];
                type_changes[i] = type_changes[j];
                type_changes[j] = tmp;
            }
        }
    }

    /*
        Assign changes to elements contained in type_groups
    */
    for (unsigned tg=0;tg<type_groups.size();tg++){
        int typegroup_id = type_groups[tg]->id;

        std::vector<int>    elements_trial;
        std::vector<int>    elements_removed;
        std::vector<int>    elements_added;
        std::vector<int>    elements_remaining;

        int monomer_id;
        for (unsigned tc=0;tc<type_changes.size();tc++) {

            monomer_id = type_changes[tc][TYPE_CHANGE_MONOMER_ID];
            if (full_type_map[monomer_id] == typegroup_id) {
                elements_removed.push_back(monomer_id);
            }

            if (type_changes[tc][TYPE_CHANGE_CHANGE_TO] == typegroup_id) {
                elements_added.push_back(monomer_id);
            }
        }

        unsigned num_added   = elements_added.size();
        unsigned num_removed = elements_removed.size();

        /*
            Set elements_remaining
        */
        elements_remaining = type_groups[tg]->elements;
        for (unsigned rm=0;rm<elements_removed.size();rm++) {
            elements_remaining.erase(std::remove(elements_remaining.begin(), elements_remaining.end(), elements_removed[rm]), elements_remaining.end());
        }

        #ifdef UNBOUND_DEBUG
        /*
            Verify the correct ordering of elements_remaining
        */
        for (unsigned elrm=0;elrm<elements_remaining.size()-1;elrm++) {
            if (elements_remaining[elrm]>= elements_remaining[elrm+1]) {
                std::cout << "Error: Unbound::eval_type_changes(): elements_remaining is uncorrectly ordered." << std::endl;
                std::exit(0);
            }
        }
        #endif

        /*
            Set elements_trial
        */
        elements_trial = elements_remaining;
        for (unsigned ad=0;ad<elements_added.size();ad++) {
            elements_trial.push_back(elements_added[ad]);
        }
        std::sort(elements_trial.begin(), elements_trial.end());

        type_groups[tg]->elements_trial     = elements_trial;
        type_groups[tg]->elements_removed   = elements_removed;
        type_groups[tg]->elements_added     = elements_added;
        type_groups[tg]->elements_remaining = elements_remaining;
    }
}



/*---------------------------------------------------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------------------------------------------------*/


double Unbound::Eval_Delta_Energy(  const arma::mat * pos_new,
                                    const arma::mat * pos_old,
                                    const std::vector<arma::ivec> & moved )
{
    double Delta_E = 0;
    int pair_move_count = gen_type_groups_moved_intervals(moved);
//    if (pair_move_count > 0)  {
    if (true) {
        for (unsigned p=0;p<pairs.size();p++) {
            #ifdef UNBOUND_DEBUG
            double dE1 = pairs[p]->Eval_Delta_Energy(pos_new,pos_old);

            double Enew = pairs[p]->Eval_Energy(pos_new);
            double Eold = pairs[p]->Eval_Energy(pos_old);

            double dE2 = Enew - Eold;
            if (dE2 > PAIR_INFINITY/10) {
                dE2 = PAIR_INFINITY;
            }
//            std::cout << "======================" << std::endl;
//            std::cout << "E_old = " << Eold << std::endl;
//            std::cout << "E_new = " << Enew << std::endl;
//
//            std::cout << "dE1 = " << dE1 << std::endl;
//            std::cout << "dE2 = " << dE2 << std::endl;

//            if (std::abs(dE1-dE2) > EQUAL_EPS) {
//                std::cout << "Error in Pair: Eval_Delta_Energy and Eval_Energy are inconsistent" << std::endl;
//                std::cout << "Part Calculation = " << dE1 << std::endl;
//                std::cout << "Full Calculation = " << dE2 << std::endl;
//                std::cout << "E_old = " << Eold << std::endl;
//                std::cout << "E_new = " << Enew << std::endl;
//
//                std::cout << arma::accu(*pos_new-*pos_old) << std::endl;
//                std::exit(0);
//            }
            Delta_E += dE2;
            #else
            Delta_E += pairs[p]->Eval_Delta_Energy(pos_new,pos_old);
            #endif
        }
    }
    return Delta_E;
}



int Unbound::gen_type_groups_moved_intervals(const std::vector<arma::ivec> & moved_intervals) {

    #ifdef UNBOUND_DEBUG
    // verify that moved_intervals is ordered
    for (unsigned i=0;i<moved_intervals.size()-1;i++){
        if (moved_intervals[i+1](EV_FROM) < moved_intervals[i](EV_FROM) && moved_intervals[i](EV_TYPE) != -1 && moved_intervals[i+1](EV_TYPE) != -1) {
            std::cout << "Error: Unbound::gen_type_groups_moved_intervals - moved_intervals is not ordered!" << std::endl;
            for (unsigned mi=0;mi<moved_intervals.size();mi++) {
                std::cout << moved_intervals[mi].t() << std::endl;
            }
            std::exit(0);
        }
    }
    #endif

//    std::cout << "####################################################################################" << std::endl;
//    std::cout << "####################################################################################" << std::endl;
//    for (unsigned mi=0;mi<moved_intervals.size();mi++) {
//        std::cout << moved_intervals[mi].t() << std::endl;
//    }

    int pair_move_count = 0;

    for (unsigned i=0;i<type_groups.size();i++){
        std::vector<int> unmoved;
        std::vector<int> moved;
        std::vector<int> individual;
        int type_group_move_count = 0;

//        std::cout << "######" << std::endl;
//        std::cout << "type group " << i << std::endl;

        for (unsigned m=0;m<moved_intervals.size();m++){
            int rge_from = moved_intervals[m](EV_FROM);
            int rge_to   = moved_intervals[m](EV_TO);
            int movtype  = moved_intervals[m](EV_TYPE);

            /*
                Discard empty placeholder intervals
            */
            if (movtype == -1) continue;

            /*
                find elements within this range
            */
            int first = -1;
            int last;

            // find first
            for (unsigned elem=0;elem<type_groups[i]->elements.size();elem++) {
                if ( rge_from <= type_groups[i]->elements[elem] ) {
                    if ( type_groups[i]->elements[elem] <= rge_to ) {
                        first = elem;
                    }
                    break;
                }
            }

            // find last
            if (first > -1) {
                last = type_groups[i]->elements.size()-1;
                for (unsigned elem=first+1;elem<type_groups[i]->elements.size();elem++) {
                    if (type_groups[i]->elements[elem] > rge_to) {
                        last = elem -1;
                        break;
                    }
                }
            }
            else {
                /*
                    No elements fall within this range
                */
                first =  0;
                last  = -1;
            }

//            std::cout << first << " - " << last << std::endl;


            if (within_EV_unmoved(movtype)) {
                unmoved.push_back(first);
                unmoved.push_back(last);
            }
            else {
                moved.push_back(first);
                moved.push_back(last);

                if (type_groups[i]->involved_in_pairinteraction) {
                    pair_move_count += last+1-first;
                    type_group_move_count++;
                }

                if (within_EV_individual(movtype)) {
                    individual.push_back(first);
                    individual.push_back(last);
                }
            }
        }
        type_groups[i]->unmoved     = unmoved;
        type_groups[i]->moved       = moved;
        type_groups[i]->individual  = individual;
        type_groups[i]->num_moved   = type_group_move_count;

//        print_vector(type_groups[i]->elements);
//        print_vector(type_groups[i]->unmoved);
//        print_vector(type_groups[i]->moved);
//        print_vector(type_groups[i]->individual);
//        std::cout << type_group_move_count << std::endl;

    }
    return pair_move_count;
}




/*---------------------------------------------------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------------------------------------------------*/



int Unbound::atom_type_exists(int id) {
    for (unsigned j=0;j<type_groups.size();j++) {
        if (type_groups[j]->id == id) {
            return true;
        }
    }
    return false;
}

int Unbound::get_index_in_type_groups(int id) {
    for (unsigned j=0;j<type_groups.size();j++) {
        if (type_groups[j]->id == id) {
            return j;
        }
    }
    return -1;
}


void Unbound::assign_type_group_elements() {
    /*
        clear current element lists
    */
    for (unsigned tg=0;tg<type_groups.size();tg++){
        type_groups[tg]->elements.clear();
        type_groups[tg]->elements.shrink_to_fit();
    }
    for (unsigned i=0;i<full_type_map.size();i++) {
//        int type         = full_type_map[i];
//        int group_index  = get_index_in_type_groups(type);
//        std::cout << i << " " << type << " " << group_index << std::endl;
        type_groups[get_index_in_type_groups(full_type_map[i])]->elements.push_back(i);
    }
}


/*---------------------------------------------------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------------------------------------------------*/
/*
                            INIT FUNCTIONS
*/
/*---------------------------------------------------------------------------------------------------------------------*/


void Unbound::read_inputfile(const std::string & pairinteraction_file) {

    InputRead * interactions_read = new InputRead(pairinteraction_file);

    // Init AtomTypes
    num_types = interactions_read->get_single_num_instances(TYPE_IDENTIFIER);
    std::vector<std::vector<std::string>> atom_str_vecs;
    for (unsigned i=0;i<num_types;i++) {
        atom_str_vecs.push_back(interactions_read->get_single_entries(i,TYPE_IDENTIFIER));
    }
    init_atoms(atom_str_vecs);

    int new_default_type = interactions_read->get_single_val<int>(DEFAULT_TYPE_IDENTIFIER);
    if (new_default_type != default_type) {
        default_type = new_default_type;
        for (unsigned i=0;i<full_type_map.size();i++) {
            full_type_map[i] = default_type;
        }
    }


    // Init Pair Interactions
    num_pairs = interactions_read->get_single_num_instances(PAIR_IDENTIFIER);
    std::vector<std::vector<std::string>> pair_str_vecs;
    for (unsigned i=0;i<num_pairs;i++) {
        pair_str_vecs.push_back(interactions_read->get_single_entries(i,PAIR_IDENTIFIER));
    }
    init_pairs(pair_str_vecs);
 }

void Unbound::init_atoms(const std::vector<std::vector<std::string>> & atom_str_vecs) {
    std::cout << std::endl << "Initializing Atom Types .." << std::endl;

    for (auto tg : type_groups)
    {
        delete tg;
    }
    type_groups.clear();
    type_groups.shrink_to_fit();
    for (unsigned i=0;i<atom_str_vecs.size();i++) {
        std::vector<double> input = strVec2doubleVec(atom_str_vecs[i]);

        // check if sufficient arguments are given
        if (input.size() < 2) {
            std::cout << " Insufficient arguments for atom type definition. Atom type definitions require at least 2 arguments: " << std::endl;
            std::cout << "  - (1) atom id " << std::endl;
            std::cout << "  - (2) hardcore repulsion radius -> set to zero for no hard core repulsion." << std::endl;
            std::exit(0);
        }
        int type_id = input[0];

        if (fmod(input[0],1)!=0) {
            std::cout << " Warning: Atom type ids need to be integer. " << input[0] << " was given. Value was set to " << type_id << std::endl;
        }

        // check if atom with given id is already initialized
        if (atom_type_exists(type_id)) {
            std::cout << " Error: atom type " << type_id << " was specified multiple times." << std::endl;
            std::exit(0);
        }
        if (type_id < 0) {
            std::cout << " Error: negative atom types are not permitted. Atom type id " << type_id << " was given." << std::endl;
            std::exit(0);
        }

        TypeGroup * newtype = new TypeGroup;

        newtype->id      = type_id;
        newtype->type.id = type_id;

        if (input[1] < 0) {
            newtype->type.hard_repulsion_cuttoff = 0;
            std::cout << " Error: negative repulsion radius specified." << std::endl;
            std::exit(0);
        }
        else {
            newtype->type.hard_repulsion_cuttoff = input[1];
        }
        std::vector<double>::const_iterator first = input.begin() + 2;
        std::vector<double>::const_iterator last  = input.end();
        std::vector<double> newVec(first, last);
        newtype->type.params = newVec;
        type_groups.push_back(newtype);
    }
    std::cout << " " << type_groups.size() << " atom types initialized" << std::endl;
}



void Unbound::init_monomer_type_assignments_from_file(const std::string & typeset_file) {

    InputRead * input = new InputRead(typeset_file);

    std::vector<std::vector<int>> range_sets;
    std::vector<std::vector<int>> individual_sets;

    int num_range = input->get_single_num_instances(INIT_TYPE_RANGE_SET_IDENFIER);
    for (unsigned i=0;i<num_range;i++) {
        range_sets.push_back(input->get_single_intvec(i,INIT_TYPE_RANGE_SET_IDENFIER));
    }

    int num_indiv = input->get_single_num_instances(INIT_TYPE_SET_IDENFIER);
    for (unsigned i=0;i<num_indiv;i++) {
        individual_sets.push_back(input->get_single_intvec(i,INIT_TYPE_SET_IDENFIER));
    }

        init_monomer_type_assignments(range_sets,individual_sets);
}

void Unbound::init_monomer_type_assignments(const std::vector<std::vector<int>> & range_sets,
                                            const std::vector<std::vector<int>> & individual_sets) {
    int type_id;
    int type_groups_id;
    int rge_from, rge_to;
    int mon_id;
    /*
        Assign ranges
    */
    for (unsigned i=0;i<range_sets.size();i++) {
        if (range_sets[i].size() < 3) {
            std::cout << " Error: Range sets for atom type init require 3 arguments (" << range_sets[i].size() << " given)." << std::endl;
            std::exit(0);
        }

        type_id         = range_sets[i][0];
        rge_from        = range_sets[i][1];
        rge_to          = range_sets[i][2];
        type_groups_id  = get_index_in_type_groups(type_id);

        if (type_groups_id < 0) {
            std::cout << " Error: attempting to assign monomer to non-existing type (" << type_id << ")." << std::endl;
            std::cout << "  Initialized types are:" << std::endl;
            for (unsigned tg=0;tg<type_groups.size();tg++) {
                std::cout << "   - " << type_groups[tg]->id << std::endl;
            }
            std::exit(0);
        }

        if (rge_from < 0) {
            std::cout << " Error: Lower bound for range sets for atom type init need to be larger or equal to 0 (" << rge_from << " given)." << std::endl;
            std::exit(0);
        }

        if (rge_to >= chain->get_num_bp()) {
            std::cout << " Error: Upper bound for range sets for atom type init exceeds number of monomers (id "
                        << rge_to << " given for " << chain->get_num_bp() << " contained monomers)." << std::endl;
            std::exit(0);
        }

        for (int rge=rge_from;rge<=rge_to;rge++) {
            full_type_map[rge] = type_id;
        }
    }

    /*
        Assign individual
    */
    for (unsigned i=0;i<individual_sets.size();i++) {
        if (individual_sets[i].size() < 2) {
            std::cout << " Error: Individual sets for atom type init require 2 arguments (" << individual_sets[i].size() << " given)." << std::endl;
            std::exit(0);
        }

        type_id         = individual_sets[i][0];
        type_groups_id  = get_index_in_type_groups(type_id);
        if (type_groups_id < 0) {
            std::cout << " Error: attempting to assign monomer to non-existing type (" << type_id << ")." << std::endl;
            std::cout << "  Initialized types are:" << std::endl;
            for (unsigned tg=0;tg<type_groups.size();tg++) {
                std::cout << "   - " << type_groups[tg]->id << std::endl;
            }
            std::exit(0);
        }
        for (unsigned j=1;j<individual_sets[i].size();j++) {
            mon_id = individual_sets[i][j];
            if (mon_id < 0 || mon_id >= chain->get_num_bp()) {
                std::cout << " Error: Invalid monomer id given for individual monomer type assignment. "
                            << "Monomer ids need to be within the range [" << 0 << "," << chain->get_num_bp()-1 << "], but "
                            << mon_id << " was given." << std::endl;
                std::exit(0);
            }

            full_type_map[mon_id] = type_id;
        }
    }
    assign_type_group_elements();
    set_backup();
}



/*---------------------------------------------------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------------------------------------------------*/













