#include "Unbound.h"


void Unbound::init_pairs(const std::vector<std::vector<std::string>> & pair_str_vecs) {
    std::cout << std::endl << "Initializing Pair Interactions .." << std::endl;

    for (auto prs : pairs)
    {
        delete prs;
    }
    pairs.clear();
    pairs.shrink_to_fit();

    for (unsigned i=0;i<pair_str_vecs.size();i++) {
        std::vector<std::string> pair_str_vec = pair_str_vecs[i];

        // check if sufficient arguments are given
        if (pair_str_vec.size() < 4) {
            std::cout << " Insufficient arguments given for pair interaction. Generic Pair interactions require at least 3 arguments" << std::endl;
            std::cout << "  - (1)  String idenfier for used potential" << std::endl;
            std::cout << "  - (2)  id for first atom type" << std::endl;
            std::cout << "  - (3)  id for second atom type" << std::endl;
            std::cout << "  - (4)  interaction cutoff" << std::endl;
            std::cout << "  - (5-) specific interaction parameters" << std::endl;
            std::exit(0);
        }

        std::string interaction_idenfier = pair_str_vec[0];
        int         atom_id1             = std::stoi(pair_str_vec[1]);
        int         atom_id2             = std::stoi(pair_str_vec[2]);
        double      cutoff               = std::stod(pair_str_vec[3]);

        std::vector<std::string> sub(&pair_str_vec[4],&pair_str_vec[pair_str_vec.size()]);
        std::vector<double> params       = strVec2doubleVec(sub);

        // check validity of arguments
        if (!atom_type_exists(atom_id1)) {
            std::cout << " Error: atom type " << atom_id1 << " specified for pair intaction '" << interaction_idenfier << "'not initialized." << std::endl;
            std::exit(0);
        }
        if (!atom_type_exists(atom_id2)) {
            std::cout << " Error: atom type " << atom_id2 << " specified for pair intaction '" << interaction_idenfier << "' not initialized." << std::endl;
            std::exit(0);
        }
        if (cutoff < 0) {
            std::cout << " Error: Negative interaction cutoff specified for pair style '" << interaction_idenfier << "'." << std::endl;
            std::exit(0);
        }

        /*---------------------------------------------------------------------------------------------------------------------*/
        /*---------------------------------------------------------------------------------------------------------------------*/
        /*---------------------------------------------------------------------------------------------------------------------*/
        /*
            New Pair Styles have to be initialized here!
        */

        if (interaction_idenfier == PAIR_TYPE_IDENFIER_LENNARD_JONES) {
            Pair * new_lj = new Pair_LJ(    chain,
                                            type_groups[get_index_in_type_groups(atom_id1)],
                                            type_groups[get_index_in_type_groups(atom_id2)],
                                            cutoff,
                                            params);
            pairs.push_back(new_lj);
        }
        else {
            std::cout << " Error: Pair Style '" << interaction_idenfier << "' unknown." << std::endl;
            std::exit(0);
        }

        /*---------------------------------------------------------------------------------------------------------------------*/
        /*---------------------------------------------------------------------------------------------------------------------*/
        /*---------------------------------------------------------------------------------------------------------------------*/

    }

    std::cout << " " << pairs.size() << " atom pair interactions initialized\n" << std::endl;
}
