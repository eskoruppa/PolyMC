#include "MultiLine.h"

bool MultiLineInstance::elem_within_bounds(unsigned elem) {
    if (elem >= singlelines.size()) {
        std::cout << "Error MultiLineInstance: For identifier '" << identifier << "'. Element " << elem << " of " << singlelines.size() << " requested" << std::endl;
        std::exit(0);
    }
    return true;
}

unsigned MultiLineInstance::amount_singlelines() {
    return singlelines.size();
}

SingleLineInstance * MultiLineInstance::get_singleline(unsigned id)  {
    elem_within_bounds(id);
    return & singlelines[id];
}

std::vector<SingleLineInstance> * MultiLineInstance::get_all_singlelines()  {
    return & singlelines;
}

bool MultiLineInstance::contains_singleline(const std::string & identifier) {
    for (unsigned i=0;i<singlelines.size();i++) {
        if (singlelines[i].identifier == identifier ) {
            return true;
        }
    }
    return false;
}





bool MultiLine::elem_within_bounds(unsigned elem) {
    if (elem >= instances.size()) {
        std::cout << "Error MultiLine: For identifier '" << identifier << "'. Element " << elem << " of " << instances.size() << " requested" << std::endl;
        std::exit(0);
    }
    return true;
}

unsigned MultiLine::amount_multilines() {
    return instances.size();
}

MultiLineInstance * MultiLine::get_instance(unsigned elem) {
    elem_within_bounds(elem);
    return & instances[elem];
}

std::vector<MultiLineInstance> * MultiLine::get_all_multilines() {
    return & instances;
}

std::vector<double> MultiLineInstance::get_single_vec(const std::string& identifier) {
    for (unsigned i=0;i<singlelines.size();i++) {
        if (singlelines[i].identifier == identifier ) {
            return singlelines[i].get_vec();
        }
    }
    std::cout << "No vector with identifier '" << identifier << "' found!" << std::endl;
    std::exit(0);
    return {-1};
}

std::vector<int>    MultiLineInstance::get_single_intvec(const std::string& identifier) {
    for (unsigned i=0;i<singlelines.size();i++) {
        if (singlelines[i].identifier == identifier ) {
            return singlelines[i].get_intvec();
        }
    }
    std::cout << "No vector with identifier '" << identifier << "' found!" << std::endl;
    std::exit(0);
    return {-1};
}

