#include "ExtraFuncs.h"

std::vector<double> strVec2doubleVec(const std::vector<std::string> & strvec)
{
    std::vector<double> vec;
    for (unsigned i=0;i<strvec.size();i++) {
        double value;
        std::istringstream ss(strvec[i]);
        ss >> value;
        vec.push_back(value);
    }
    return vec;
}

std::vector<int> strVec2intVec(const std::vector<std::string> & strvec)
{
    std::vector<int> vec;
    for (unsigned i=0;i<strvec.size();i++) {
        int value;
        std::istringstream ss(strvec[i]);
        ss >> value;
        vec.push_back(value);
    }
    return vec;
}
