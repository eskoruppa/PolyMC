#include "IDB.h"

/*
TODO:
- Allow for variable type definition in Database file (indicated by x)
- Allow for spaces in monomer-types definition

*/

IDB::IDB(const std::string & filename, bool all_oligomers_required)
: fn(filename), all_olis_required(all_oligomers_required) {
    read_file(fn);
}

IDB::~IDB() {
}

//////////////////////////////////////////////////////////
///////////////////// ACCESSORS //////////////////////////
//////////////////////////////////////////////////////////


double IDB::get_disc_len() {
    return disc_len;
}

unsigned IDB::get_interaction_range() {
    return interaction_range;
}

bool IDB::get_avg_inconsist() {
    return avg_inconsist;
}

bool IDB::IDBU_contained(const std::string & type) {
    for (unsigned i=0;i<idbus.size();i++) {
        if (idbus[i].type == type) {
            return true;
        }
    }
    return false;
}

// IDBU * IDB::get_IDBU(const std::string & type) {
    // for (unsigned i=0;i<idbus.size();i++) {
    //     if (idbus[i].type == type) {
    //         return &idbus[i];
    //     }
    // }
    // std::cout << "Error IDB::get_IDBU - requested type '" << type << "' was not provided in IDB file." << std::endl;
    // std::exit(0);
// }

IDBU * IDB::get_IDBU(const std::string & type) {

    // std::cout << "Requesting ologomer '" << type << "'." << std::endl;

    // check if idbu is already contained
    for (unsigned i=0;i<idbus.size();i++) {
        if (idbus[i].type == type) {
            return &idbus[i];
        }
    }

    if (!input->contains_multi(type)) {
        throw std::runtime_error("IDB::get_IDBU: missing entry.");
    }

    MultiLine * ml = input->get_multiline(type);
    if (ml->amount_multilines() > 1) {
        throw std::runtime_error("IDB::get_IDBU: Double mention of interaction in IDB file.");
    }

    MultiLineInstance * mli = ml->get_instance(0);
    unsigned num_single = mli->amount_singlelines();

    if (num_single != sets_per_unit+1) {
        throw std::runtime_error("IDB::get_IDBU: Invalid number of subentries.");
    }

    IDBU * new_idbu = new IDBU();
    // IDBU new_idbu;
    new_idbu->type=type;
    std::vector<SingleLineInstance> singles = *(mli->get_all_singlelines());
    for (unsigned j=0;j<num_single;j++) {
        if (singles[j].identifier == IDB_ID_GROUNDSTATE) {
            new_idbu->Theta0 = arma::colvec(singles[j].get_vec());
        }
        else {
            new_idbu->methods.push_back(singles[j].identifier);
            new_idbu->interactions.push_back(singles[j].get_vec());
        }
    }
    idbus.push_back(*new_idbu);
    return new_idbu;
}

//////////////////////////////////////////////////////////
///////////////////// READ IDB ///////////////////////////
//////////////////////////////////////////////////////////

bool IDB::read_file(const std::string & filename) {
    // InputRead idbinput(filename);
    InputRead * inputread = new InputRead(filename);
    input = inputread;

    bp_types          = input->get_single_val<std::string>(IDB_ID_MONOMER_TYPES);
    disc_len          = input->get_single_val<double>(IDB_ID_DISC_LEN);
    avg_inconsist     = input->get_single_val<bool>(IDB_ID_AVG_INCONSIST);

    interaction_range = input->get_single_val<unsigned>(IDB_ID_INTERACTION_RANGE);
    sets_per_unit     = 1+2*interaction_range;

    // gen_typeperm();
    // gen_idbus(&input);
    std::vector<IDBU>().swap(idbus);

    return true;
}



void IDB::gen_idbus(InputRead * input) {

    bool consistent_inputfile = true;
    unsigned missing_counter  = 0;
    std::vector<IDBU>().swap(idbus);

    /*
        acquiring core couplings

        These have to be contained otherwise an Error will be raised and the program terminates here.
    */
    for (unsigned i=0;i<typeperms_core.size();i++) {
        if (!input->contains_multi(typeperms_core[i])) {
            if (all_olis_required) {
                std::cout << "Error: Missing interaction in IDB file: " << typeperms_core[i] << std::endl;
                consistent_inputfile=false;
            }
            missing_counter++;
            continue;
        }

        MultiLine * ml = input->get_multiline(typeperms_core[i]);
        if (ml->amount_multilines() > 1) {
            std::cout << "Error: Double mention of interaction in IDB file: " << typeperms_core[i] << std::endl;
            consistent_inputfile=false;
            continue;
        }

        MultiLineInstance * mli = ml->get_instance(0);
        unsigned num_single = mli->amount_singlelines();

        if (num_single != sets_per_unit+1) {
            std::cout << "Error: interaction entry " << typeperms_core[i] << " in IDB file should contain " << sets_per_unit+1  << " subentries, but " << num_single << " entries were found!" << std::endl;
            consistent_inputfile=false;
            continue;
        }

        IDBU new_idbu;
        new_idbu.type=typeperms_core[i];
        std::vector<SingleLineInstance> singles = *(mli->get_all_singlelines());
        for (unsigned j=0;j<num_single;j++) {
            if (singles[j].identifier == IDB_ID_GROUNDSTATE) {
                new_idbu.Theta0 = arma::colvec(singles[j].get_vec());
            }
            else {
                new_idbu.methods.push_back(singles[j].identifier);
                new_idbu.interactions.push_back(singles[j].get_vec());
            }
        }
        idbus.push_back(new_idbu);
    }
    /*
        acquiring termini couplings

        These do not necessarily have to be contained. If not contained the first core partial match (with x substituted for any type) will be assigned.
    */

    for (unsigned i=0;i<typeperms_termini.size();i++) {
        if (input->contains_multi(typeperms_termini[i])) {
            /*
                If Terminus interaction is contained in IDB file
            */

            MultiLine * ml = input->get_multiline(typeperms_termini[i]);
            if (ml->amount_multilines() > 1) {
                std::cout << "Error: Double mention of interaction in IDB file: " << typeperms_termini[i] << std::endl;
                consistent_inputfile=false;
                continue;
            }

            MultiLineInstance * mli = ml->get_instance(0);
            unsigned num_single = mli->amount_singlelines();

            if (num_single != sets_per_unit+1) {
                std::cout << "Error: interaction entry " << typeperms_termini[i] << " in IDB file should contain " << sets_per_unit+1  << " subentries, but " << num_single << " entries were found!" << std::endl;
                consistent_inputfile=false;
                continue;
            }

            IDBU new_idbu;
            new_idbu.type=typeperms_termini[i];
            std::vector<SingleLineInstance> singles = *(mli->get_all_singlelines());
            for (unsigned j=0;j<num_single;j++) {
                if (singles[j].identifier == IDB_ID_GROUNDSTATE) {
                    new_idbu.Theta0 = arma::colvec(singles[j].get_vec());
                }
                else {
                    new_idbu.methods.push_back(singles[j].identifier);
                    new_idbu.interactions.push_back(singles[j].get_vec());
                }
            }
            idbus.push_back(new_idbu);
        }
        else {
            /*
                If Terminus interaction is not contained in IDB file
            */
            std::string placeholder;
            for (unsigned j=0;j<typeperms_core.size();j++) {
                if (type_terminus_match(typeperms_core[j],typeperms_termini[i])) {
                    placeholder = typeperms_core[j];
                    break;
                }
            }
            if (!input->contains_multi(placeholder)) {
//                std::cout << "Error: Missing interaction in IDB file: " << placeholder << std::endl;
                consistent_inputfile=false;
                continue;
            }

            MultiLine * ml = input->get_multiline(placeholder);
            MultiLineInstance * mli = ml->get_instance(0);
            unsigned num_single = mli->amount_singlelines();

            IDBU new_idbu;
            new_idbu.type=typeperms_termini[i];
            std::vector<SingleLineInstance> singles = *(mli->get_all_singlelines());
            for (unsigned j=0;j<num_single;j++) {
                if (singles[j].identifier == IDB_ID_GROUNDSTATE) {
                    new_idbu.Theta0 = arma::colvec(singles[j].get_vec());
                }
                else {
                    new_idbu.methods.push_back(singles[j].identifier);
                    new_idbu.interactions.push_back(singles[j].get_vec());
                }
            }
            idbus.push_back(new_idbu);
        }
    }

    if (!consistent_inputfile) {
        std::cout << "Program terminated because of inconsistent IDB file" << std::endl;
        std::cout << missing_counter << " couplings are missing!" << std::endl;
        std::exit(0);
    }
}



void IDB::gen_typeperm() {
    long long int num_perms;
    if (bp_types.length()>1) {
        num_perms = std::pow(bp_types.length(),(interaction_range+1)*2) + 2*bp_types.length()*(std::pow(bp_types.length(),(interaction_range+1)*2-1) - 1)/(bp_types.length()-1);
    }
    else {
        num_perms = 1+2*((interaction_range+1)*2-1);
    }
    if (num_perms>IDB_MAX_PERMUATIONS) {
        std::cout << "Error IDB::gentypeperm - possible number of interaction permuations exceeds " << IDB_MAX_PERMUATIONS << " ("  << num_perms << ")" << std::endl;
        std::exit(0);
    }

    std::vector<std::string>().swap(typeperms);
    std::vector<std::string>().swap(typeperms_core);
    std::vector<std::string>().swap(typeperms_termini);

    typepermtree("",(interaction_range+1)*2,&typeperms_core);

    std::string seq = "";
    // for (unsigned i=1;i<(interaction_range+1)*2;i++) {
    for (unsigned i=1;i<=interaction_range;i++) {
        seq+="x";
        typepermtree(seq,(interaction_range+1)*2,&typeperms_termini);
    }
    unsigned terminussize = typeperms_termini.size();
    for (unsigned i=0;i<terminussize;i++) {
        if (typeperms_termini[i][0]=='x') {
            std::string rev = typeperms_termini[i];
            std::reverse(rev.begin(), rev.end());
            typeperms_termini.push_back(rev);
        }
    }

    typeperms = typeperms_core;
    typeperms.insert( typeperms.end(), typeperms_termini.begin(), typeperms_termini.end() );

//    for (int i=0;i<typeperms.size();i++) {
//        std::cout << typeperms[i] << std::endl;
//    }
//    std::cout << typeperms_core.size() << std::endl;
//    std::cout << typeperms_termini.size() << std::endl;
//
//    std::cout << typeperms.size() << std::endl;
//    std::cout << num_perms << std::endl;
//    std::exit(0);
}

void IDB::typepermtree(std::string seq, unsigned total,std::vector<std::string> * permlist) {
    if (seq.length() == total) permlist->push_back(seq);
    else {
        for (unsigned i=0;i<bp_types.length();i++) {
            typepermtree(seq+bp_types[i],total,permlist);
        }
    }
}

bool IDB::type_terminus_match(std::string A, std::string B) {
    if (A.size() != B.size()) {
        return false;
    }
    for (unsigned i=0;i<A.size();i++) {
        if (A[i] != B[i] && A[i] != 'x' && B[i] != 'x') {
            return false;
        }
    }
    return true;
}



