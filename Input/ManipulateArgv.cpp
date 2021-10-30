#include "ManipulateArgv.h"
#include <stdexcept>
#include <iostream>



std::vector<std::string> add2argv(const std::vector<std::string> & argv, const std::string & flag, const std::vector<std::string> & flag_args) {
    if (flag[0] != '-') {
        throw std::invalid_argument("flag nees to start with '-'");
    }
    std::vector<std::string> nargv = argv;
    // remove flag if already contained
    for (int ii=0;ii<nargv.size();ii++) {
        if (nargv[ii] == flag) {
            int first = ii;
            int last  = nargv.size();
            for (int jj=ii+1;jj<nargv.size();jj++) {
                std::string str = nargv[jj];
                if (nargv[jj].at(0) == '-') {
                    last = jj;
                    break;
                }
            }
            nargv.erase(nargv.begin()+first,nargv.begin()+last);
            break;
        }
    }
    nargv.push_back(flag);
    for (int ii=0;ii<flag_args.size();ii++) {
        nargv.push_back(flag_args[ii]);
    }
    return nargv;
}

std::vector<std::string> add2argv(const std::vector<std::string> & argv, const std::string & flag, const std::string & flag_arg) {
    if (flag.at(0) != '-') {
        throw std::invalid_argument("flag nees to start with '-'");
    }
    std::vector<std::string> nargv = argv;
    // remove flag if already contained
    for (int ii=0;ii<nargv.size();ii++) {
        if (nargv[ii] == flag) {
            int first = ii;
            int last  = nargv.size();
            for (int jj=ii+1;jj<nargv.size();jj++) {
                std::string str = nargv[jj];
                if (nargv[jj].at(0) == '-') {
                    last = jj;
                    break;
                }
            }
            nargv.erase(nargv.begin()+first,nargv.begin()+last);
            break;
        }
    }
    nargv.push_back(flag);
    nargv.push_back(flag_arg);
    return nargv;
}
