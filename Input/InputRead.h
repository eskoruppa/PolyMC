#ifndef __INPUT_READ_INCLUDED__
#define __INPUT_READ_INCLUDED__

#include <string>    // string
#include <iostream>  // input output
#include <stdexcept> // exceptions
#include <fstream>   // read file
#include <sstream>
#include <vector>
#include <armadillo> // armadillo
#include <algorithm>

#include "SingleLine.h"
#include "MultiLine.h"


#define IS_MULTI     true
#define IS_NOT_MULTI false



class InputRead;

class InputRead {

protected:
    std::string filename;
    std::vector<SingleLine> singles;
    std::vector<MultiLine>  multies;
    std::vector<std::vector<std::string>> linelists;
    bool infilefound;


public:

    InputRead(const std::string& filename);
    ~InputRead();

////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////// CALL METHODS ////////////////////////////////////////////

public:

    bool filefound();
    void testfilefound();

/////////////////////////////// SINGLES ///////////////////////////////

    bool     contains_single(const std::string& identifier);
    unsigned get_single_num_instances(const std::string& identifier);
    std::vector<std::string> get_single_idlist();

    SingleLine * get_singleline(const std::string& identifier);

    template<typename T>
    T get_single_val(const std::string& identifier);
    template<typename T>
    T get_single_val(unsigned elem, const std::string& identifier);

    template<typename T>
    T get_single_val(const std::string& identifier,T default_val);
    template<typename T>
    T get_single_val(unsigned elem,const std::string& identifier,T default_val);

    std::vector<double> get_single_vec(const std::string& identifier);
    std::vector<double> get_single_vec(unsigned elem, const std::string& identifier);

    std::vector<int>    get_single_intvec(const std::string& identifier);
    std::vector<int>    get_single_intvec(unsigned elem, const std::string& identifier);

    std::vector<std::string> get_single_entries(const std::string& identifier);
    std::vector<std::string> get_single_entries(unsigned elem, const std::string& identifier);
    std::vector<std::string> get_single_full(const std::string& identifier);
    std::vector<std::string> get_single_full(unsigned elem, const std::string& identifier);

 /////////////////////////////// MULTIES ///////////////////////////////

    bool     contains_multi(const std::string& identifier);
    unsigned get_multi_num_instances(const std::string& identifier);
    std::vector<std::string> get_multi_idlist();

    MultiLine * get_multiline(const std::string& identifier);

//    bool     multi_contains_subentry(unsigned elem, const std::string& identifier,const std::string& subentry);
//    unsigned multi_amount_subentries(unsigned elem, const std::string& identifier);
//    unsigned multi_amount_specific_subentries(unsigned elem, const std::string& identifier,const std::string& subentry);
//
//    template<typename T>
//    T get_multi_entry(unsigned elem,const std::string& identifier, const std::string& entry_identifier, T default_val,unsigned entry_elem);
//
////    template<typename T>
//    std::vector<double> get_multi_entry_vec(unsigned elem,const std::string& identifier, const std::string& entry_identifier ,unsigned entry_elem);
//    std::vector<std::string> get_multi_entry_strvec(unsigned elem,const std::string& identifier, const std::string& entry_identifier ,unsigned entry_elem);
//
//    std::vector<std::string> get_multi_full_entry(unsigned elem,const std::string& identifier, unsigned entry_num);
//    std::vector<std::string> get_multi_full_entry(const std::string& identifier, unsigned entry_num);


////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////// READ FILE //////////////////////////////////////////////

protected:
    void read_file();

    void init_entries(const std::vector<std::vector<std::string>> &linelists);

    bool is_multientry(unsigned id, const std::vector<std::vector<std::string>> &linelists);
    bool is_subentry(unsigned id, const std::vector<std::vector<std::string>> &linelists);
    bool entry_already_initiated(const std::string &identifier, bool multi);
    void find_all_entries(const std::string &identifier, const std::vector<std::vector<std::string>> &linelists, bool multi);

    void print_entries();

    void get_linelist(std::string &str, std::vector<std::string> &output);
    void split_string(std::string &str, const std::string &delims, std::vector<std::string> &output);

};

////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////// TEMPLATES //////////////////////////////////////////////


template<typename T>
T InputRead::get_single_val(const std::string& identifier) {
    return get_single_val<T>(0,identifier);
}

template<typename T>
T InputRead::get_single_val(unsigned elem, const std::string& identifier) {
    for (unsigned i=0;i<singles.size();i++) {
        if (identifier == singles[i].identifier) {
            return singles[i].get_first<T>(elem);
        }
    }
    std::cout << " Error identifier '" << identifier << "' not found in input file '" << filename << "'" << std::endl;
    std::exit(0);
}

template<typename T>
T InputRead::get_single_val(const std::string& identifier,T default_val) {
    return get_single_val<T>(0,identifier,default_val);
}

template<typename T>
T InputRead::get_single_val(unsigned elem,const std::string& identifier,T default_val) {
    for (unsigned i=0;i<singles.size();i++) {
        if (identifier == singles[i].identifier) {
            return singles[i].get_first<T>(elem);
        }
    }
    return default_val;
}


//template<typename T>
//T InputRead::get_multi_entry(unsigned elem, const std::string& identifier, const std::string& entry_identifier, T default_val,unsigned entry_elem) {
//    for (unsigned i=0;i<multi_entries.size();i++) {
//        if (identifier == multi_entries[i].identifier) {
//            if (multi_entries[i].entries.size() <= elem) {
//                std::cout << " Error: multi entry with identifier '" << identifier << "' contains only " << single_entries[i].entries.size() << " elements" << std::endl;
//                std::cout << "          Element " << elem << " was requested!" << std::endl;
//                std::exit(0);
//            }
//
//            unsigned counter = 0;
//            for (unsigned j=0;j<multi_entries[i].entries[elem].size();j++) {
//                if (entry_identifier == multi_entries[i].entries[elem][j].identifier) {
//                    if (counter==entry_elem) {
//
//                        T value;
//                        std::istringstream ss(multi_entries[i].entries[elem][j].entries[0]);
//                        ss >> value;
//                        return value;
//
//                    }
//                    counter++;
//                }
//            }
//
//            std::cout << " Error: sub_entry '" << entry_identifier << "' of multi entry with identifier '" << identifier << "' contains only " << counter << " elements" << std::endl;
//            std::cout << "          Element " << entry_elem << " was requested!" << std::endl;
//            std::exit(0);
//        }
//    }
//    return default_val;
//}


#endif
