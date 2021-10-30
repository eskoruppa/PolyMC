#include "SingleLine.h"


//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////
/////////// SingleLineInstance ///////////////////////////////

unsigned SingleLineInstance::get_num_entries() {
    return entries.size();
}

std::vector<double> SingleLineInstance::get_vec(){
    std::vector<double> vec;
    for (unsigned i=0;i<entries.size();i++) {
        double value;
        std::istringstream ss(entries[i]);
        ss >> value;
        vec.push_back(value);
    }
    return vec;
}

std::vector<int> SingleLineInstance::get_intvec(){
    std::vector<int> vec;
    for (unsigned i=0;i<entries.size();i++) {
        int value;
        std::istringstream ss(entries[i]);
        ss >> value;
        vec.push_back(value);
    }
    return vec;
}


std::vector<std::string> SingleLineInstance::get_entries(){
    return entries;
}

std::vector<std::string> SingleLineInstance::get_full(){
    std::vector<std::string> vec;
    vec.push_back(identifier);
    for (unsigned i=0;i<entries.size();i++) {
        vec.push_back(entries[i]);
    }
    return vec;
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////
/////////// SingleLine ///////////////////////////////////////


bool SingleLine::elem_within_bounds(unsigned elem) {
    if (elem >= instances.size()) {
        std::cout << "Error SingleLine: For identifier '" << identifier << "'. Element " << elem << " of " << instances.size() << " requested" << std::endl;
        std::exit(0);
    }
    return true;
}

std::vector<double> SingleLine::get_vec(unsigned elem) {
    elem_within_bounds(elem);
    return instances[elem].get_vec();
}

std::vector<double> SingleLine::get_vec() {
    return get_vec(0);
}

std::vector<std::string> SingleLine::get_entries(unsigned elem) {
    elem_within_bounds(elem);
    return instances[elem].get_entries();
}

std::vector<std::string> SingleLine::get_entries() {
    return get_entries(0);
}

std::vector<std::string> SingleLine::get_full(unsigned elem) {
    elem_within_bounds(elem);
    return instances[elem].get_full();
}

std::vector<std::string> SingleLine::get_full() {
    return get_full(0);
}















