#ifndef __MULTILINE__
#define __MULTILINE__

#include <string>    // string
#include <iostream>  // input output
#include <stdexcept> // exceptions
#include <fstream>   // read file
#include <sstream>
#include <vector>
#include <armadillo> // armadillo
#include <algorithm>

#include "SingleLine.h"


struct MultiLineInstance;
struct MultiLineInstance {
    std::string identifier;
    std::vector<SingleLineInstance> singlelines;

    bool elem_within_bounds(unsigned elem);

    unsigned amount_singlelines();
    SingleLineInstance * get_singleline(unsigned elem);
    std::vector<SingleLineInstance> * get_all_singlelines();

    bool contains_singleline(const std::string & identifier);

    template<typename T>
    T get_single_val(const std::string& identifier);
    template<typename T>
    T get_single_val(const std::string& identifier, T default_val);

    std::vector<double> get_single_vec(const std::string& identifier);
//    std::vector<double> get_single_vec(unsigned elem, const std::string& identifier);
//
    std::vector<int>    get_single_intvec(const std::string& identifier);
//    std::vector<int>    get_single_intvec(unsigned elem, const std::string& identifier);

};


template<typename T>
T MultiLineInstance::get_single_val(const std::string& identifier) {
    for (unsigned i=0;i<singlelines.size();i++) {
        if (singlelines[i].identifier == identifier ) {
            return singlelines[i].get_first<T>();
        }
    }
    std::cout << "Identifier '" << identifier << "' not found!" << std::endl;
    std::exit(0);
    return -1;
}

template<typename T>
T MultiLineInstance::get_single_val(const std::string& identifier, T default_val) {
    for (unsigned i=0;i<singlelines.size();i++) {
        if (singlelines[i].identifier == identifier ) {
            return singlelines[i].get_first<T>();
        }
    }
    return default_val;
}



struct MultiLine;
struct MultiLine {
    std::string identifier;
    std::vector<MultiLineInstance> instances;

    bool elem_within_bounds(unsigned elem);

    unsigned amount_multilines();
    MultiLineInstance * get_instance(unsigned elem);
    std::vector<MultiLineInstance> * get_all_multilines();
};


#endif

