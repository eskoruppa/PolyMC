#ifndef __GENINPUTFILE_INCLUDED__
#define __GENINPUTFILE_INCLUDED__

#include <string>    // string
#include <iostream>  // input output
#include <stdexcept> // exceptions
#include <fstream>   // read file
#include <sstream>
#include <vector>
#include <armadillo> // armadillo
#include <algorithm>


#define GENINFILE_SIMSETUP      "Simulation Setup"
#define GENINFILE_SEEDS         "Seed"
#define GENINFILE_INTERACTIONS  "Interactions"
#define GENINFILE_PARAMS        "Simulation Parameters"
#define GENINFILE_OUTPUT        "Output"
#define GENINFILE_DUMPS         "Dumps"

#define GENINPUTFILE_SEPERATOR_LENGTH 60
#define GENINPUTFILE_LEFTSEP 2


class GenInputFile;
struct InputFileCategory;
struct InpitFileEntry;



class GenInputFile {

protected:
    std::vector<InputFileCategory*>   categories;


public:
    GenInputFile();
    ~GenInputFile();

    void add_category(const std::string & name);
    void add_category(const std::string & name,const std::string & explanation);
    void add_category_explantion(const std::string & category_name,const std::string & explanation);

    bool add_entry(const std::string & category_name, const std::string & identifier, const std::string & value);
    template<typename T>
    bool add_entry(const std::string & category_name, const std::string & identifier, T value);
    void add_newline(const std::string & category_name);

    void generate_input_file(const std::string & filename);

};

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

struct InputFileCategory {

public:
    std::string name;
    std::string explanation;
    std::vector<std::vector<std::string>> entries;

public:
    InputFileCategory(const std::string & name);
    InputFileCategory(const std::string & name, const std::string & explanation);
    ~InputFileCategory();

    void set_explanation(const std::string & explanation);
    bool add_entry(const std::string & identifier, const std::string & value);

    std::vector<std::string> get_lines();

};

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

template<typename T>
bool GenInputFile::add_entry(const std::string & category_name, const std::string & identifier, T value) {
    for (unsigned i=0;i<categories.size();i++) {
        if (categories[i]->name == category_name) {
            std::ostringstream ss;
            ss << value;
            return categories[i]->add_entry(identifier,ss.str());
        }
    }
    std::cout << "template" << std::endl;
    std::string ecptstr = "Inputfile category '" + category_name + "' does not exists!";
    throw std::invalid_argument( ecptstr );
}




#endif
