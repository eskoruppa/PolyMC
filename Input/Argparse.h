#ifndef __INPUT_ARGPARSER_INCLUDED__
#define __INPUT_ARGPARSER_INCLUDED__

#include <string>    // string
//#include <stdexcept> // exceptions
//#include <fstream>   // read file
#include <sstream>
#include <vector>


template<typename T>
T convstr(T var, const std::string &argstr) {
    T value;
    std::istringstream ss(argstr);
    ss >> value;
    return value;
}

template<typename T>
T parse_arg(T default_val,const std::string &keyword, const std::vector<std::string> & argv) {
    int argc = argv.size();
    std::string argstr;
    for (int i=1;i<argc;i++) {
		argstr = argv[i];
		if (argstr==keyword && argc > i+1) {
            argstr = argv[i+1];
            T value;
            std::istringstream ss(argstr);
            ss >> value;
            return value;
		}
    }
    return default_val;
}

std::string parse_arg(const std::string& default_val,const std::string &keyword, const std::vector<std::string> & argv);
bool        parse_flag(const std::string &keyword, const std::vector<std::string> & argv);
bool        parse_flag(bool default_val,const std::string &keyword, const std::vector<std::string> & argv);
bool        parse_flag_atpos(const std::string &keyword, int flag_pos, const std::vector<std::string> & argv);
std::string parse_str_atpos(int pos, const std::vector<std::string> & argv);


//template<typename T>
//T parse_arg(T default_val,const std::string &keyword, int argc, const char ** argv) {
//    std::string argstr;
//    for (int i=1;i<argc;i++) {
//		argstr = argv[i];
//		if (argstr==keyword && argc > i+1) {
//            argstr = argv[i+1];
//            T value;
//            std::istringstream ss(argstr);
//            ss >> value;
//            return value;
//		}
//    }
//    return default_val;
//}
//
//std::string parse_arg(const std::string& default_val,const std::string &keyword, int argc, const char ** argv);
//bool        parse_flag(const std::string &keyword, int argc, const char ** argv);
//bool        parse_flag(bool default_val,const std::string &keyword, int argc, const char ** argv);
//bool        parse_flag_atpos(const std::string &keyword, int flag_pos , int argc, const char ** argv);
//std::string parse_str_atpos(int pos , int argc, const char ** argv);

#endif
