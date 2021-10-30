#include "ReadFunctions.h"

void get_linelist(std::string &str, std::vector<std::string> &output) {
    output.clear();
    if (str[0] == ' ' || str[0] == '\t'){
        output.push_back("\t");
        str.erase(0, 1);
    }
    std::replace( str.begin(), str.end(), '\t', ' ');
    while (str[0] == ' '){
        str.erase(0, 1);
    }
    std::vector<std::string> tmp;
    std::string delims{ " =" };
    split_string(str,delims,tmp);

    output.insert(output.end(),tmp.begin(),tmp.end());
}

void split_string(std::string &str, const std::string &delims, std::vector<std::string> &output)
{
    size_t beg, pos = 0;
    while ((beg = str.find_first_not_of(delims, pos)) != std::string::npos)
    {
        pos = str.find_first_of(delims, beg + 1);
        output.push_back(str.substr(beg, pos - beg));
    }
}


