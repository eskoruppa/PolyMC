#ifndef __INTERACTIONDATABASE_INCLUDED__
#define __INTERACTIONDATABASE_INCLUDED__

#include <armadillo>

#include <iostream> // std::cout
#include <string>   // std::string
#include <algorithm>
#include <stdexcept>
#include <fstream>  // std::ifstream
#include <vector>
using namespace std;


struct IDBUnit;
class IDB;


struct IDBUnit {
std::string type;
arma::colvec Theta0;
arma::cube   Ms;
};


class IDB {

private:
    std::string     fn;
    IDBUnit *       units;
    unsigned int    amount_units;
    unsigned int    sets_per_unit;
    unsigned int    interaction_range;
    std::string     bp_types;
    double          disc_len;

public:

    IDB(std::string filename);
    ~IDB();
    bool read_file();

    double   get_disc_len();
    unsigned get_interaction_range();
    bool     IDBU_contained(const std::string& type, int till=-1);
    IDBUnit* IDBU_get_unit(const std::string& type);
    bool     IDBU_get_terminus(IDBUnit& term_unit,const std::string& vartype);

protected:
    bool IDB_vartype_match(const std::string& vartype, const std::string& type);

private:
    bool IDBU_add_unit(const std::string& type, const arma::colvec& Theta0, const arma::cube& Ms,int id);

    std::string       line_get_first_word(const std::string& line);
    std::string       line_get_arg (const std::string& line,unsigned int num);
    std::string       line_merge_tail(const std::string& line,unsigned int num);
    arma::mat    line_get_mat (const std::string& line,unsigned int start);
    arma::colvec line_get_vec (const std::string& line,unsigned int start);

    void split_string(const std::string &str, std::vector<std::string> &output);
};

#endif
