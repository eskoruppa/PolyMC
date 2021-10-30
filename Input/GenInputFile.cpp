#include "GenInputFile.h"


GenInputFile::GenInputFile()
{

}

GenInputFile::~GenInputFile() {
    for (unsigned i=0;i<categories.size();i++) {
        delete categories[i];
    }
}

void GenInputFile::add_category(const std::string & name) {
    // check if category already exists
    for (unsigned i=0;i<categories.size();i++) {
        if (categories[i]->name == name) {
            throw std::invalid_argument( "Inputfile category already exists!" );
        }
    }
    InputFileCategory * category = new InputFileCategory(name);
    categories.push_back(category);
}

void GenInputFile::add_category(const std::string & name,const std::string & explanation) {
    // check if category already exists
    for (unsigned i=0;i<categories.size();i++) {
        if (categories[i]->name == name) {
            throw std::invalid_argument( "Inputfile category already exists!" );
        }
    }
    InputFileCategory * category = new InputFileCategory(name,explanation);
    categories.push_back(category);
}

void GenInputFile::add_category_explantion(const std::string & category_name,const std::string & explanation) {
    for (unsigned i=0;i<categories.size();i++) {
        if (categories[i]->name == category_name) {
            categories[i]->set_explanation(explanation);
            return;
        }
    }
    std::string ecptstr = "Inputfile category '" + category_name + "' does not exists!";
    throw std::invalid_argument( ecptstr );
}

bool GenInputFile::add_entry(const std::string & category_name, const std::string & identifier, const std::string & value) {
    for (unsigned i=0;i<categories.size();i++) {
        if (categories[i]->name == category_name) {
            return categories[i]->add_entry(identifier,value);
        }
    }
    std::string ecptstr = "Inputfile category '" + category_name + "' does not exists!";
    throw std::invalid_argument( ecptstr );
}

void GenInputFile::add_newline(const std::string & category_name) {
    for (unsigned i=0;i<categories.size();i++) {
        if (categories[i]->name == category_name) {
            categories[i]->add_entry("","");
            return;
        }
    }
    std::string ecptstr = "Inputfile category '" + category_name + "' does not exists!";
    throw std::invalid_argument( ecptstr );
}


void GenInputFile::generate_input_file(const std::string & filename) {
    int sep = GENINPUTFILE_SEPERATOR_LENGTH;
    int add = GENINPUTFILE_LEFTSEP*2 + 2;
    for (unsigned i=0;i<categories.size();i++) {
        int entrysep = categories[i]->name.length() + add;
        if (entrysep > sep) {
            sep = entrysep;
        }
    }

    std::string sepstr(sep,'#');
    sepstr += "\n";
    std::string leftsep(GENINPUTFILE_LEFTSEP,'#');

    std::ofstream infile;
    infile.open(filename, std::ofstream::out | std::ofstream::trunc);
    for (unsigned i=0;i<categories.size();i++) {
        std::string rightsep(sep - categories[i]->name.length()-2-GENINPUTFILE_LEFTSEP,'#');
        infile << sepstr;
        infile << leftsep << " " << categories[i]->name << " " << rightsep << "\n";
        infile << sepstr;
//      TODO: SPLIT WORDS CORRECTLY IN EXPLANATION
//        if (categories[i]->explanation.length() > 0) {
//            int explen = sep - 3;
//            for (unsigned l=0;l<categories[i]->explanation.length()/explen+1;l++) {
//                int left = categories[i]->explanation.length() - l*explen;
//                if (left > explen) {
//                    left = explen;
//                }
//                infile << "# " << categories[i]->explanation.substr(l*explen, left) << "\n";
//            }
//        }
        infile << "\n";
        std::vector<std::string> lines = categories[i]->get_lines();
        for (unsigned j=0;j<lines.size();j++) {
            infile << lines[j] << "\n";
        }
        infile << "\n";

    }
    infile.close();
}



////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


InputFileCategory::InputFileCategory(const std::string & name)
{
    this->name = name;
    this->explanation = "";

}
InputFileCategory::InputFileCategory(const std::string & name, const std::string & explanation)
{
    this->name = name;
    this->explanation = explanation;
}

InputFileCategory::~InputFileCategory() {}

void InputFileCategory::set_explanation(const std::string & explanation) {
    this->explanation = explanation;
}

bool InputFileCategory::add_entry(const std::string & identifier, const std::string & value) {
    if (identifier != "") {
        for (unsigned i=0;i<entries.size();i++) {
            if (entries[i][0] == identifier) {
                entries[i][1] = value;
                return true;
            }
        }
    }
    entries.push_back({identifier,value});
    return false;
}

std::vector<std::string> InputFileCategory::get_lines() {
    int longest_identifier = 0;
    for (unsigned i=0;i<entries.size();i++) {
        if (entries[i][0].length() > longest_identifier) {
            longest_identifier = entries[i][0].length();
        }
    }
    std::vector<std::string> lines;
    for (unsigned i=0;i<entries.size();i++) {
        int len_idenfifier = entries[i][0].length();
        if (len_idenfifier==0) {
            lines.push_back("");
        }
        else {
            std::string connect(longest_identifier-len_idenfifier, ' ');
            connect += " = ";
            lines.push_back(entries[i][0]+connect+entries[i][1]);
        }
    }
    return lines;
}

