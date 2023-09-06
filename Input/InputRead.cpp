#include "InputRead.h"


InputRead::InputRead(const std::string& filename)
: filename(filename)
{
    read_file();
    init_entries(linelists);
//    print_entries();
}

InputRead::~InputRead() {}


////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////// CALL METHODS ////////////////////////////////////////////


bool InputRead::filefound() {
    return infilefound;
}

void InputRead::testfilefound() {
    if (!infilefound) {
        std::cout << "Error: Input file not found!" << std::endl;
        std::exit(0);
    }
}


////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////// SINGLE METHODS ///////////////////////////////////////////


bool InputRead::contains_single(const std::string& identifier) {
    for (unsigned i=0;i<singles.size();i++) {
        if (singles[i].identifier == identifier) {
            return true;
        }
    }
    return false;
}

unsigned InputRead::get_single_num_instances(const std::string& identifier) {
    for (unsigned i=0;i<singles.size();i++) {
        if (identifier == singles[i].identifier) {
            return singles[i].instances.size();
        }
    }
    return 0;
}

std::vector<std::string> InputRead::get_single_idlist() {
    std::vector<std::string> elems;
    for (unsigned i=0;i<singles.size();i++) {
        elems.push_back(singles[i].identifier);
    }
    return elems;
}

SingleLine * InputRead::get_singleline(const std::string& identifier) {
    for (unsigned i=0;i<singles.size();i++) {
        if (identifier == singles[i].identifier) {
            return &singles[i];
        }
    }
    std::cout << "Error singleline '" << identifier << "' not found." << std::endl;
    std::exit(0);
}

std::vector<double> InputRead::get_single_vec(const std::string& identifier) {
    return get_single_vec(0,identifier);
}

std::vector<double> InputRead::get_single_vec(unsigned elem, const std::string& identifier) {
    for (unsigned i=0;i<singles.size();i++) {
        if (identifier == singles[i].identifier) {
            return singles[i].get_vec(elem);
        }
    }
    std::cout << "Error: InputRead::get_single_vec couldn't find element with identifier '" << identifier << "'." << std::endl;
    std::exit(0);
}

std::vector<int> InputRead::get_single_intvec(const std::string& identifier) {
    return get_single_intvec(0,identifier);
}

std::vector<int> InputRead::get_single_intvec(unsigned elem, const std::string& identifier) {
    for (unsigned i=0;i<singles.size();i++) {
        if (identifier == singles[i].identifier) {
            std::vector<double> DoubleVec = singles[i].get_vec(elem);
            std::vector<int> intVec(DoubleVec.begin(), DoubleVec.end());
            return intVec;
        }
    }
    std::cout << "Error: InputRead::get_single_vec couldn't find element with identifier '" << identifier << "'." << std::endl;
    std::exit(0);
}


std::vector<std::string> InputRead::get_single_entries(const std::string& identifier) {
    return get_single_entries(0,identifier);
}

std::vector<std::string> InputRead::get_single_entries(unsigned elem, const std::string& identifier) {
    for (unsigned i=0;i<singles.size();i++) {
        if (identifier == singles[i].identifier) {
            return singles[i].get_entries(elem);
        }
    }
    std::cout << "Error: InputRead::get_single_entries couldn't find element with identifier '" << identifier << "'." << std::endl;
    std::exit(0);
}

std::vector<std::string> InputRead::get_single_full(const std::string& identifier) {
    return get_single_full(0,identifier);
}

std::vector<std::string> InputRead::get_single_full(unsigned elem, const std::string& identifier) {
    for (unsigned i=0;i<singles.size();i++) {
        if (identifier == singles[i].identifier) {
            return singles[i].get_full(elem);
        }
    }
    std::cout << "Error: InputRead::get_single_full couldn't find element with identifier '" << identifier << "'." << std::endl;
    std::exit(0);
}



////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////// MUTLI METHODS ///////////////////////////////////////////


bool InputRead::contains_multi(const std::string& identifier) {
    for (unsigned i=0;i<multies.size();i++) {
        if (multies[i].identifier == identifier) {
            return true;
        }
    }
    return false;
}

unsigned InputRead::get_multi_num_instances(const std::string& identifier) {
    for (unsigned i=0;i<multies.size();i++) {
        if (identifier == multies[i].identifier) {
            return multies[i].instances.size();
        }
    }
    return 0;
}

std::vector<std::string> InputRead::get_multi_idlist() {
    std::vector<std::string> elems;
    for (unsigned i=0;i<multies.size();i++) {
        elems.push_back(multies[i].identifier);
    }
    return elems;
}

MultiLine * InputRead::get_multiline(const std::string& identifier) {
    for (unsigned i=0;i<multies.size();i++) {
        if (identifier == multies[i].identifier) {
            return &multies[i];
        }
    }
    std::cout << "Error: InputRead::get_multiline couldn't find element with identifier '" << identifier << "'." << std::endl;
    std::exit(0);
}


//bool InputRead::multi_contains_subentry(unsigned elem, const std::string& identifier,const std::string& subentry) {
//    for (unsigned i=0;i<multi_entries.size();i++) {
//        if (identifier == multi_entries[i].identifier) {
//            for (unsigned j=0;j<multi_entries[i].entries[elem].size();j++) {
//                if (multi_entries[i].entries[elem][j].identifier == subentry) {
//                    return true;
//                }
//            }
//        }
//    }
//    return false;
//}
//
//unsigned InputRead::multi_amount_subentries(unsigned elem, const std::string& identifier){
//    for (unsigned i=0;i<multi_entries.size();i++) {
//        if (identifier == multi_entries[i].identifier) {
//            return multi_entries[i].entries[elem].size();
//        }
//    }
//    return 0;
//}
//
//unsigned InputRead::multi_amount_specific_subentries(unsigned elem, const std::string& identifier,const std::string& subentry) {
//    for (unsigned i=0;i<multi_entries.size();i++) {
//        if (identifier == multi_entries[i].identifier) {
//            unsigned counter = 0;
//            for (unsigned j=0;j<multi_entries[i].entries[elem].size();j++) {
//                if (multi_entries[i].entries[elem][j].identifier == subentry) {
//                    counter++;
//                }
//            }
//            return counter;
//        }
//    }
//    return 0;
//}
//
//std::vector<double> InputRead::get_multi_entry_vec(unsigned elem,const std::string& identifier, const std::string& entry_identifier ,unsigned entry_elem) {
//    for (unsigned i=0;i<multi_entries.size();i++) {
//        if (identifier == multi_entries[i].identifier) {
//            if (multi_entries[i].entries.size() <= elem) {
//                std::cout << " Error: multi entry with identifier '" << identifier << "' contains only " << single_entries[i].entries.size() << " elements" << std::endl;
//                std::cout << "          Element " << elem << " was requested!" << std::endl;
//                std::exit(0);
//            }
//            unsigned counter = 0;
//            for (unsigned j=0;j<multi_entries[i].entries[elem].size();j++) {
//                if (entry_identifier == multi_entries[i].entries[elem][j].identifier) {
//                    if (counter==entry_elem) {
//
//                        std::vector<double> vec;
//                        for (unsigned k=0;k<multi_entries[i].entries[elem][j].entries.size();k++) {
//                            double value;
//                            std::istringstream ss(multi_entries[i].entries[elem][j].entries[k]);
//                            ss >> value;
//                            vec.push_back(value);
//                        }
//                        return vec;
//                    }
//                    counter++;
//                }
//            }
//            std::cout << " Error: sub_entry '" << entry_identifier << "' of multi entry with identifier '" << identifier << "' contains only " << counter << " elements" << std::endl;
//            std::cout << "          Element " << entry_elem << " was requested!" << std::endl;
//            std::exit(0);
//        }
//    }
//    std::cout << "Error: InputRead::get_multi_entry_vec couldn't find requested entry" << std::endl;
//    std::exit(0);
//}
//
//std::vector<std::string> InputRead::get_multi_entry_strvec(unsigned elem,const std::string& identifier, const std::string& entry_identifier ,unsigned entry_elem) {
//    for (unsigned i=0;i<multi_entries.size();i++) {
//        if (identifier == multi_entries[i].identifier) {
//            if (multi_entries[i].entries.size() <= elem) {
//                std::cout << " Error: multi entry with identifier '" << identifier << "' contains only " << single_entries[i].entries.size() << " elements" << std::endl;
//                std::cout << "          Element " << elem << " was requested!" << std::endl;
//                std::exit(0);
//            }
//            unsigned counter = 0;
//            for (unsigned j=0;j<multi_entries[i].entries[elem].size();j++) {
//                if (entry_identifier == multi_entries[i].entries[elem][j].identifier) {
//                    if (counter==entry_elem) {
//
//                        std::vector<std::string> vec;
//                        for (unsigned k=0;k<multi_entries[i].entries[elem][j].entries.size();k++) {
//                            vec.push_back(multi_entries[i].entries[elem][j].entries[k]);
//                        }
//                        return vec;
//                    }
//                    counter++;
//                }
//            }
//            std::cout << " Error: sub_entry '" << entry_identifier << "' of multi entry with identifier '" << identifier << "' contains only " << counter << " elements" << std::endl;
//            std::cout << "          Element " << entry_elem << " was requested!" << std::endl;
//            std::exit(0);
//        }
//    }
//    std::cout << "Error: InputRead::get_multi_entry_vec couldn't find requested entry" << std::endl;
//    std::exit(0);
//}
//
//
//std::vector<std::string> InputRead::get_multi_full_entry(unsigned elem,const std::string& identifier, unsigned entry_num) {
//    for (unsigned i=0;i<multi_entries.size();i++) {
//        if (identifier == multi_entries[i].identifier) {
//            if (multi_entries[i].entries.size() <= elem) {
//                std::cout << " Error: multi entry with identifier '" << identifier << "' contains only " << single_entries[i].entries.size() << " elements" << std::endl;
//                std::cout << "          Element " << elem << " was requested!" << std::endl;
//                std::exit(0);
//            }
//
//            if (multi_entries[i].entries[elem].size() <= entry_num) {
//                std::cout << " Error: Requested entry number " << entry_num << " of multi entry with identifier '" << identifier << "', which only contains " << multi_entries[i].entries[elem].size() << " entries" << std::endl;
//                std::exit(0);
//            }
//
//            std::vector<std::string> vec;
//            vec.push_back(multi_entries[i].entries[elem][entry_num].identifier);
//            for (unsigned k=0;k<multi_entries[i].entries[elem][entry_num].entries.size();k++) {
//                vec.push_back(multi_entries[i].entries[elem][entry_num].entries[k]);
//            }
//            return vec;
//        }
//    }
//    std::cout << "Error: InputRead::get_multi_entry_vec couldn't find requested entry" << std::endl;
//    std::exit(0);
//
//}
//
//std::vector<std::string> InputRead::get_multi_full_entry(const std::string& identifier, unsigned entry_num) {
//    return get_multi_full_entry(0,identifier,entry_num);
//}


////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////// READ FILE //////////////////////////////////////////////

void InputRead::read_file() {

    std::ifstream file;
    std::string   line;

    file.open(filename);
    if (file.is_open())
    {
        while(!file.eof()) {
            std::getline(file,line);
            std::vector<std::string> linelist;
            get_linelist(line,linelist);
            if (linelist.size() > 0) {
                if (linelist.size() == 1 && linelist[0] == "\t") {
                    continue;
                }
                std::string first;
                for (unsigned i=0;i<linelist.size();i++) {
                    if (linelist[i] != "\t" && linelist[i] != " ") {
                        first = linelist[i][0];
                        break;
                    }
                }
                if (first!="#") {
                    linelists.push_back(linelist);
                }
            }
        }
        infilefound=true;
    }
    else {
        infilefound=false;
        std::cout << "\nERROR:" << std::endl;
        std::cout << " Unable to find specified input file '" << filename << "'" << std::endl;
        std::cout << "simulation terminated .. " << std::endl;
        std::exit(0);
    }
    file.close();
}

////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////



void InputRead::init_entries(const std::vector<std::vector<std::string>> &linelists) {

    /*
        Init Single Line Entries
    */
    std::string identifier;
    for (unsigned i=0;i<linelists.size();i++) {
        if (!is_subentry(i,linelists)) {
            identifier = linelists[i][0];
            find_all_entries(identifier, linelists, is_multientry(i, linelists));
        }
    }
}


bool InputRead::is_multientry(unsigned id, const std::vector<std::vector<std::string>> &linelists) {
    if (id < linelists.size()-1 && linelists[id+1][0] == "\t") {
        return true;
    }
    return false;
}

bool InputRead::is_subentry(unsigned id, const std::vector<std::vector<std::string>> &linelists) {
    if (linelists[id][0] == "\t")
        return true;
    return false;
}

bool InputRead::entry_already_initiated(const std::string &identifier, bool multi) {
    if (multi) {
        for (unsigned i=0;i<multies.size();i++) {
            if (identifier==multies[i].identifier)
                return true;
        }
    }
    else {
        for (unsigned i=0;i<singles.size();i++) {
            if (identifier==singles[i].identifier)
                return true;
        }
    }
    return false;
}

void InputRead::find_all_entries(const std::string &identifier, const std::vector<std::vector<std::string>> &linelists, bool multi) {
    if (multi) {
        if (!entry_already_initiated(identifier,IS_MULTI)) {

            std::vector<MultiLineInstance> mls;
            for (unsigned i=0;i<linelists.size();i++) {
                if (!is_subentry(i, linelists) && is_multientry(i, linelists) && linelists[i][0] == identifier) {

                    std::vector<SingleLineInstance> singlelines;
                    for (unsigned j=i+1;j<linelists.size();j++) {
                        if (is_subentry(j,linelists)) {

                            std::string sl_identifier = linelists[j][1];
                            std::vector<std::string>::const_iterator first = linelists[j].begin() + 2;
                            std::vector<std::string>::const_iterator last  = linelists[j].end();
                            std::vector<std::string> newVec(first, last);

                            SingleLineInstance sl;
                            sl.identifier = sl_identifier;
                            sl.entries    = newVec;
                            singlelines.push_back(sl);
                        }
                        else {
                            break;
                        }
                    }

                    MultiLineInstance newml;
                    newml.identifier  = identifier;
                    newml.singlelines = singlelines;
                    mls.push_back(newml);
                }
            }

            MultiLine ml;
            ml.identifier = identifier;
            ml.instances = mls;
            multies.push_back(ml);
        }
    }
    else {
        if (!entry_already_initiated(identifier,IS_NOT_MULTI)) {
            std::vector<SingleLineInstance> singlelines;
            for (unsigned i=0;i<linelists.size();i++) {
                if (!is_subentry(i, linelists) && !is_multientry(i, linelists) && linelists[i][0] == identifier) {
                    std::vector<std::string>::const_iterator first = linelists[i].begin() + 1;
                    std::vector<std::string>::const_iterator last  = linelists[i].end();
                    std::vector<std::string> newVec(first, last);

                    SingleLineInstance sl;
                    sl.identifier = identifier;
                    sl.entries    = newVec;
                    singlelines.push_back(sl);
                }
            }
            SingleLine SL;
            SL.identifier  = identifier;
            SL.instances = singlelines;
            singles.push_back(SL);
        }
    }
}


void InputRead::print_entries() {
    std::cout << "##################################" << std::endl;
    std::cout << "----------SINGLE-ENTRIES----------" << std::endl;
    std::cout << "##################################" << std::endl;

    for (unsigned i=0;i<singles.size();i++) {
        std::cout << "ID: " << singles[i].identifier << std::endl;
        for (unsigned j=0;j<singles[i].instances.size();j++) {
            std::cout << "--------" << std::endl;
            std::cout << "entry " << j << ":" << std::endl;
            for (unsigned k=0;k<singles[i].instances[j].entries.size();k++) {
                std::cout << singles[i].instances[j].entries[k] << std::endl;
            }
        }
        std::cout << "##################################" << std::endl;
    }

    std::cout << "num multi entries = " << multies.size() << std::endl;

    std::cout << "##################################" << std::endl;
    std::cout << "----------MULTI-ENTRIES-----------" << std::endl;
    std::cout << "##################################" << std::endl;

    for (unsigned i=0;i<multies.size();i++) {
        std::cout << "ID: " << multies[i].identifier << std::endl;
        std::cout << " " << multies[i].instances.size() << " instances found" << std::endl;
        for (unsigned j=0;j<multies[i].instances.size();j++) {
            std::cout << "---------------" << std::endl;
            std::cout << "entry " << j << ":" << std::endl;

            for (unsigned k=0;k<multies[i].instances[j].singlelines.size();k++) {
                std::cout << "subID: " << multies[i].instances[j].singlelines[k].identifier << std::endl;
                for (unsigned l=0;l<multies[i].instances[j].singlelines[k].entries.size();l++) {
                    std::cout << " " << multies[i].instances[j].singlelines[k].entries[l] << std::endl;
                }
            }
        }
        std::cout << "##################################" << std::endl;
    }
}


////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////// SPLIT LINE IN LINELIST /////////////////////////////////////////


void InputRead::get_linelist(std::string &str, std::vector<std::string> &output) {
    
    int hashid = str.find("#");
    if (hashid != std::string::npos) {
        if (hashid == 0) {
            str = "";
        }
        else {
            str = str.substr (0,hashid); 
        }
    }
    
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

void InputRead::split_string(std::string &str, const std::string &delims, std::vector<std::string> &output)
{
    size_t beg, pos = 0;
    while ((beg = str.find_first_not_of(delims, pos)) != std::string::npos)
    {
        pos = str.find_first_of(delims, beg + 1);
        output.push_back(str.substr(beg, pos - beg));
    }
}

////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////



/*
TEST InputRead


int main(int argc, const char **argv) {

    string inputfn = "testin";
    inputfn  = parse_arg(inputfn, "-in",argc,argv);
    InputRead input(inputfn);
    input.testfilefound();

//
    cout << "#################################" << endl;
    cout << "#################################" << endl;

    cout << input.contains_single("T") << endl;
    cout << input.get_single_num_instances("T") << endl;

    std::vector<std::string> single_idlist = input.get_single_idlist();
    for (unsigned i=0;i<single_idlist.size();i++){
        cout << single_idlist[i] << endl;
    }

    double def = 3.1415;
    cout << input.get_single_val(4,"T",def) << endl;
    std::string defs = "hello";
    cout << input.get_single_val(0,"str",defs) << endl;

    cout << input.get_single_val(2,"str",defs) << endl;
    cout << input.get_single_val(3,"str",defs) << endl;
    std::string id = "vec";

    cout << endl;

    vector<double> vec1,vec2;

    vec1 = input.get_single_vec("vec");
    for (int i=0;i<vec1.size();i++){
        cout << vec1[i] << endl;
    }

    cout << endl;
    vec2 = input.get_single_vec(2,"vec");
    for (int i=0;i<vec1.size();i++){
        cout << vec2[i] << endl;
    }

    vector<string> vec3,vec4;
    cout << endl;
    vec3 = input.get_single_entries("vec");
    for (int i=0;i<vec3.size();i++){
        cout << vec3[i] << endl;
    }

    cout << endl;
    vec4 = input.get_single_entries(2,"vec");
    for (int i=0;i<vec4.size();i++){
        cout << vec4[i] << endl;
    }

    cout << endl;
    vec3 = input.get_single_full("vec");
    for (int i=0;i<vec3.size();i++){
        cout << vec3[i] << endl;
    }

    cout << endl;
    vec4 = input.get_single_full(2,"vec");
    for (int i=0;i<vec4.size();i++){
        cout << vec4[i] << endl;
    }


    cout << "#################################" << endl;
    cout << "#################################" << endl;

    cout << input.contains_multi("PIV") << endl;
    cout << input.get_multi_num_instances("PIV") << endl;

    vector<string> idlist = input.get_multi_idlist();
    for (int i=0;i<idlist.size();i++){
        cout << idlist[i] << endl;
    }

    MultiLine * ml = input.get_multiline("PIV");

    unsigned num = ml->amount_multilines();
    cout << num << "instances of " << ml->identifier << "found" << endl;

    MultiLineInstance * mli = ml->get_instance(2);
    std::vector<SingleLineInstance> * allsls = mli->get_all_singlelines();

    for (unsigned j=0;j<(*allsls).size();j++){
        vector<string> v = (*allsls)[j].get_full();
        cout << "~~~~~~~~~~~~~" << endl;
        for (unsigned k=0;k<v.size();k++) {
            cout << v[k] << endl;
        }
    }




    std::vector<MultiLineInstance> * mlis = ml->get_all_multilines();



//    cout << input.multi_contains_subentry(0,"PIV","F") << endl;
//    cout << input.multi_contains_subentry(1,"PIV","E") << endl;
//    cout << input.multi_contains_subentry(2,"PIV","F") << endl;
//
//    cout << "#############" << endl;
//
//    cout << input.multi_amount_subentries(0,"PIV") << endl;
//    cout << input.multi_amount_subentries(1,"PIV") << endl;
//    cout << input.multi_amount_subentries(2,"PIV") << endl;
//    cout << input.multi_amount_subentries(0,"ROT") << endl;
//
//    cout << "#############" << endl;
//
//    cout << input.multi_amount_specific_subentries(0,"PIV","F") << endl;
//    cout << input.multi_amount_specific_subentries(1,"PIV","E") << endl;
//    cout << input.multi_amount_specific_subentries(2,"PIV","F") << endl;
//
//
//    cout << "#############" << endl;
//
//    cout << input.get_multi_entry(1,"PIV", "E", 2.5 ,1) << endl;
//    cout << input.get_multi_entry(1,"PIV", "E", 2.5 ,2) << endl;
//    cout << input.get_multi_entry(1,"PIV", "E", 2.5 ,0) << endl;
//
//    cout << "#############" << endl;
//    cout << "#############" << endl;
//    cout << "#############" << endl;
//
//    vector<double> vec3 = input.get_multi_entry_vec(0,"ROT","H234asd",0);
//    cout << "len vec3 = " << vec3.size() << endl;
//    for (int i=0;i<vec3.size();i++){
//        cout << vec3[i] << endl;
//    }
//
//    cout << "#############" << endl;
//
//    vector<double> vec4 = input.get_multi_entry_vec(0,"ROT","H234asd",1);
//    cout << "len vec4 = " << vec4.size() << endl;
//    for (int i=0;i<vec4.size();i++){
//        cout << vec4[i] << endl;
//    }
//
//    cout << "#############" << endl;
//
//    vector<double> vec5 = input.get_multi_entry_vec(0,"ROT","H234asd",2);
//    cout << "len vec5 = " << vec5.size() << endl;
//    for (int i=0;i<vec5.size();i++){
//        cout << vec5[i] << endl;
//    }


}
*/




