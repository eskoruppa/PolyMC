#include "../PolyMC.h"

std::string PolyMC::read_seq(unsigned num_bp) {

    if (seq_fn == "") {
        std::cout << "\nWarning: No sequence file specified. The sequence was intialized as poly 'a' .." << std::endl;
        seq = "";
        for (unsigned i=0;i<num_bp;i++) {
            seq += "a";
        }
        return seq;
    }
    else {
        std::ifstream seqfile;
        std::string   line;

        seqfile.open(seq_fn);
        if (!seqfile.is_open()) {
            std::cout << "\nError: Provided sequence file could not be opened" << std::endl;
            std::exit(0);
        }

        std::vector<std::string> first_words;
        std::string first;
        while(!seqfile.eof()) {
            std::getline(seqfile,line);
            first = get_first_word(line);
            if (first != "") {
                first_words.push_back(first);
            }
        }
        seqfile.close();

        if (first_words.size() == 0) {
            std::cout << "\nError: Provided sequence file is emtpy" << std::endl;
            std::exit(0);
            return "";
        }

        if (first_words.size() == 1) {
            seq = first_words[0];
            if (seq.length() != num_bp) {
                std::cout << "\nError: Provided sequence length (" << seq.length() << ") does not match specified number of monomers (" << num_bp << ")." << std::endl;
                std::exit(0);
            }
            return seq;
        }

        std::string seqmode = first_words[0];
        std::string partial = first_words[1];
        unsigned partlen = partial.length();

        if (seqmode=="all" || seqmode=="All" || seqmode=="match" || seqmode=="Match" || seqmode=="exact" || seqmode=="Exact") {
            if (partlen != num_bp) {
                std::cout << "\nError: Provided sequence length (" << partial.length() << ") does not match specified number of monomers (" << num_bp << ")." << std::endl;
                std::exit(0);
            }
            return partial;

        }
        if (seqmode=="repeat" || seqmode=="Repeat" || seqmode=="poly" || seqmode=="Poly") {

            if (partlen == num_bp) {
                return partial;
            }

            if (partlen > num_bp) {
                std::cout << "\nError: Provided partial sequence is longer than number of monomers" << std::endl;
                std::exit(0);
            }

            seq = "";
            unsigned num = num_bp/partlen;
            for (unsigned i=0;i<num;i++) {
                seq += partial;
            }
            unsigned rest = num_bp%partlen;
            if (rest > 0) {
                seq += partial.substr(0, rest);
            }
            return seq;
        }

        std::cout << "\nError: Provided sequence mode " << seqmode << " is unknown." << std::endl;
        std::exit(0);
    }
}

std::string PolyMC::get_first_word(std::string & str) {
    std::vector<std::string> str_list = get_str_list(str);
    if (str_list.size()>0) {
        return str_list[0];
    }
    return "";
}

std::vector<std::string> PolyMC::get_str_list(std::string &str) {
    std::replace( str.begin(), str.end(), '\t', ' ');
    while (str[0] == ' '){
        str.erase(0, 1);
    }
    std::vector<std::string> output;
    std::string delims = " ";
    split_string(str,delims,output);
    return output;
}

void PolyMC::split_string(std::string &str, const std::string &delims, std::vector<std::string> &output)
{
    size_t beg, pos = 0;
    while ((beg = str.find_first_not_of(delims, pos)) != std::string::npos)
    {
        pos = str.find_first_of(delims, beg + 1);
        output.push_back(str.substr(beg, pos - beg));
    }
}
