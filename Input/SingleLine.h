#ifndef __SINGLELINE__
#define __SINGLELINE__

#include <string>    // string
#include <iostream>  // input output
#include <stdexcept> // exceptions
#include <fstream>   // read file
#include <sstream>
#include <vector>
#include <armadillo> // armadillo
#include <algorithm>


struct SingleLineInstance;
struct SingleLineInstance {
    std::string identifier;
    std::vector<std::string> entries;

    unsigned get_num_entries();

    template<typename T>
    T get_first();
    template<typename T>
    T get_nth(unsigned nth);

    std::vector<double> get_vec();
    std::vector<int>    get_intvec();
    std::vector<std::string> get_entries();
    std::vector<std::string> get_full();
};

template<typename T>
T SingleLineInstance::get_first() {
    return get_nth<T>(0);
}

template<typename T>
T SingleLineInstance::get_nth(unsigned nth) {
    if (nth >= entries.size() ) {
        std::cout << "Error SingleLineInstance: in Instance '" << identifier << "'. Element " << nth << " of "  << entries.size() << " requested" << std::endl;
        std::exit(0);
    }
    T value;
    std::istringstream ss(entries[nth]);
    ss >> value;
    return value;
}



struct SingleLine;
struct SingleLine {
    std::string identifier;
    std::vector<SingleLineInstance> instances;
    bool elem_within_bounds(unsigned elem);

    unsigned get_num_instances();
    template<typename T>
    T get_first();
    template<typename T>
    T get_first(unsigned elem);

    template<typename T>
    T get_nth(unsigned nth);
    template<typename T>
    T get_nth(unsigned elem,unsigned nth);

    std::vector<double> get_vec(unsigned elem);
    std::vector<double> get_vec();
    std::vector<std::string> get_entries(unsigned elem);
    std::vector<std::string> get_entries();
    std::vector<std::string> get_full(unsigned elem);
    std::vector<std::string> get_full();
};


template<typename T>
T SingleLine::get_first() {
    return get_first<T>(0);
}

template<typename T>
T SingleLine::get_first(unsigned elem) {
    elem_within_bounds(elem);
    return instances[elem].get_first<T>();
}

template<typename T>
T SingleLine::get_nth(unsigned nth) {
    return get_nth<T>(0,nth);
}
template<typename T>
T SingleLine::get_nth(unsigned elem,unsigned nth) {
    elem_within_bounds(elem);
    return instances[elem].get_nth<T>(nth);
}



#endif

