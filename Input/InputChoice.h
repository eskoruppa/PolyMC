#ifndef __INPUT_CHOICE_INCLUDED__
#define __INPUT_CHOICE_INCLUDED__

#include "InputRead.h"
#include "Argparse.h"

#include <string>
#include <sstream>

template<typename T>
T InputChoice_get_single(const std::string &keyword, InputRead * input, const std::vector<std::string> argv, T default_val) {
    T val;
    if (input && input->filefound()) {
        val = input->get_single_val<T> (keyword,default_val);
        default_val = val;
    }
    std::string key = keyword;
    if (key[0] != '-') {
        key = "-"+key;
    }
    val = parse_arg(default_val,key,argv);
    return val;
}

template<typename T>
T InputChoice_get_single(const std::vector<std::string> &keywords, InputRead * input, const std::vector<std::string> argv, T default_val) {
    T val = default_val;
    for (unsigned i=0;i<keywords.size();i++) {
        if (input && input->filefound()) {
            val = input->get_single_val<T> (keywords[i],default_val);
            default_val = val;
        }
    }
    std::string key;
    for (unsigned i=0;i<keywords.size();i++) {
        key = keywords[i];
        if (key[0] != '-') {
            key = "-"+key;
        }
        val = parse_arg(val,key,argv);
    }
    return val;
}


//template<typename T>
//T InputChoice_get_single(const std::string &keyword, InputRead * input, int argc, const char ** argv, T default_val) {
//    T val;
//    if (input && input->filefound()) {
//        val = input->get_single_val<T> (keyword,default_val);
//        default_val = val;
//    }
//    std::string key = keyword;
//    if (key[0] != '-') {
//        key = "-"+key;
//    }
//    val = parse_arg(default_val,key,argc,argv);
//    return val;
//}
//
//template<typename T>
//T InputChoice_get_single(const std::vector<std::string> &keywords, InputRead * input, int argc, const char ** argv, T default_val) {
//    T val = default_val;
//    for (unsigned i=0;i<keywords.size();i++) {
//        if (input && input->filefound()) {
//            val = input->get_single_val<T> (keywords[i],default_val);
//            default_val = val;
//        }
//    }
//    std::string key;
//    for (unsigned i=0;i<keywords.size();i++) {
//        key = keywords[i];
//        if (key[0] != '-') {
//            key = "-"+key;
//        }
//        val = parse_arg(val,key,argc,argv);
//    }
//    return val;
//}

#endif
