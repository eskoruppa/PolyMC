#ifndef __IDB_INCLUDED__
#define __IDB_INCLUDED__

#include <armadillo>

#include <iostream> // std::cout
#include <string>   // std::string
#include <algorithm>
#include <stdexcept>
#include <fstream>  // std::ifstream
#include <vector>

#include "../Input/InputRead.h"

#define IDB_ID_INTERACTION_RANGE "interaction_range"
#define IDB_ID_MONOMER_TYPES     "monomer_types"
#define IDB_ID_DISC_LEN          "discretization"
#define IDB_ID_AVG_INCONSIST     "avg_inconsist"
#define IDB_ID_TYPE              "type"
#define IDB_ID_GROUNDSTATE       "vec"

#define IDB_MAX_PERMUATIONS      1e7

#define IDB_ID_GROUNDSTATE       "vec"



struct IDBU;
class IDB;

struct IDBU {
std::string                         type;
arma::colvec                        Theta0;
std::vector<std::vector<double>>    interactions;
std::vector<std::string>            methods;
};


class IDB {

private:
    std::string       fn;
    std::vector<IDBU> idbus;

    unsigned int      amount_units;
    unsigned int      sets_per_unit;

    std::vector<std::string> typeperms;
    std::vector<std::string> typeperms_core;
    std::vector<std::string> typeperms_termini;

    unsigned int                interaction_range;
    std::string                 bp_types;
    double                      disc_len;
    bool                        avg_inconsist = false;

public:

    IDB(const std::string & filename);
    ~IDB();

public:
    double   get_disc_len();
    unsigned get_interaction_range();
    bool     get_avg_inconsist();
    bool     IDBU_contained(const std::string & type);
    IDBU *   get_IDBU(const std::string & type);



private:
    bool read_file(const std::string & filename);
    void gen_idbus(InputRead * input);
    void gen_typeperm();
    void typepermtree(std::string seq, unsigned total, std::vector<std::string> * permlist);
    bool type_terminus_match(std::string A, std::string B);



//    bool     IDBU_contained(const std::string& type, int till=-1);
//    IDBU   * IDBU_get_unit(const std::string& type);
//    bool     IDBU_get_terminus(IDBU& term_unit,const std::string& vartype);
//
//protected:
//    bool IDB_vartype_match(const std::string& vartype, const std::string& type);
//
//private:
//    bool IDBU_add_unit(const std::string& type, const arma::colvec& Theta0, const arma::cube& Ms,int id);
//
//    std::string       line_get_first_word(const std::string& line);
//    std::string       line_get_arg (const std::string& line,unsigned int num);
//    std::string       line_merge_tail(const std::string& line,unsigned int num);
//    arma::mat    line_get_mat (const std::string& line,unsigned int start);
//    arma::colvec line_get_vec (const std::string& line,unsigned int start);
//
//    void split_string(const std::string &str, std::vector<std::string> &output);
};

#endif
