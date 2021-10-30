#ifndef __MANIPULATE_ARGV_INCLUDED__
#define __MANIPULATE_ARGV_INCLUDED__
#include <vector>
#include <string>


std::vector<std::string> add2argv(const std::vector<std::string> & argv, const std::string & flag, const std::vector<std::string> & flag_args);
std::vector<std::string> add2argv(const std::vector<std::string> & argv, const std::string & flag, const std::string & flag_arg);

#endif
