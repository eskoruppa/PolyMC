#include "InteractionDatabase.h"

/*
TODO:
- Allow for variable type definition in Database file (indicated by x)
- Allow for spaces in monomer-types definition

*/

IDB::IDB(string filename)
: fn(filename) {
    read_file();
}

IDB::~IDB() {
    delete[] units;
}

//////////////////////////////////////////////////////////
///////////////////// ACCESSORS //////////////////////////
//////////////////////////////////////////////////////////


double IDB::get_disc_len() {
    return disc_len;
}

unsigned IDB::get_interaction_range() {
    return interaction_range;
}

bool IDB::IDBU_contained(const string& type, int till) {
    if (till==-1) {till=amount_units;}
    for (int i=0;i<till;++i) {
        if (units[i].type==type) {
            return true;
        }
    }
    return false;
}

IDBUnit* IDB::IDBU_get_unit(const string& type) {
    for (unsigned int i=0;i<amount_units;++i) {
        if (units[i].type==type) {
            return &units[i];
        }
    }
    throw std::domain_error("IDB::IDBU_get_unit(): Trying to access non-existing interaction type");
    return nullptr;
}

bool IDB::IDBU_get_terminus(IDBUnit& term_unit, const string& vartype) {
/*
left=true means that the left terminius in the sequence is chosen. This function then looks for the
average over all partial agreements of type sequences over all found left hand continuations.
*/
    term_unit.type   = vartype;
    term_unit.Theta0 = arma::zeros(3);
    term_unit.Ms     = arma::zeros(3,3,sets_per_unit);
    int counter=0;

    for (unsigned int i=0;i<amount_units;++i) {
        if (IDB_vartype_match(vartype,units[i].type)) {
            term_unit.Theta0 = term_unit.Theta0 + units[i].Theta0;
            term_unit.Ms     = term_unit.Ms     + units[i].Ms;
            ++counter;
        }
    }
    if (counter==0) {
        cout << "Warning: IDB::IDBU_get_terminus: No vartype match found for " << vartype << ". All IDBUnit elements set to zero. " << endl;
        return false;
    }
    term_unit.Theta0 = term_unit.Theta0/counter;
    term_unit.Ms     = term_unit.Ms/counter;
    return true;
}


bool IDB::IDB_vartype_match(const string& vartype, const string& type) {
    if (vartype.length()!=type.length()) {
        return false;
    }
    for (unsigned i=0;i<type.length();++i) {
        if (vartype[i]!=type[i] && vartype[i] != 'x') {
            return false;
        }
    }
    return true;
}



//////////////////////////////////////////////////////////
//////// REATOUT FILE AND CONSTRUCT IDB UNITS ////////////
//////////////////////////////////////////////////////////

bool IDB::IDBU_add_unit(const string& type, const arma::colvec& Theta0, const arma::cube& Ms,int id) {
    if (IDBU_contained(type,id)) {
        cout << "Warning: IDB::IDBU_add_unit(): trying to add existing unit: " << type << endl;
        return false;
    }
    if (type.length()!=(interaction_range+1)*2) {
        cout << "Warning: IDB::IDBU_add_unit(): interaction '"<< type << "' is inconsistence with specified range ("<< interaction_range <<"). Interaction ignored." << endl;
        return false;
    }
    for (unsigned int i=0;i<type.length();++i) {
        bool consistent=false;
        for (unsigned int j=0;j<bp_types.length();++j) {
            if (type[i]==bp_types[j]) {
                consistent=true;
                break;
            }
        }
        if (!consistent) {
            cout << "Warning: IDB::IDBU_add_unit(): interaction '"<< type << "' contains unexpected type. Interaction ignored." << endl;
            return false;
        }
    }
    IDBUnit* unit = new IDBUnit;
    unit->type    = type;
    unit->Theta0  = Theta0;
    unit->Ms      = Ms;
    units[id] = *unit;
    delete unit;
    return true;
}

bool IDB::read_file() {

    bool rang_read = false;
    bool flav_read = false;
    bool disc_read = false;

    string rang_id = "interaction-range";
    string flav_id = "monomer-types";
    string disc_id = "discretization";

    string stiff_id  = "mat";
    string ground_id = "vec";

    string line;
    string identifier;
    ifstream file;

    cout << " - Loading Interaction Database - " << endl;
    cout << "opening file  " << fn << endl;

    file.open(fn);
    if (file.is_open())
    {
        bool reading_init = true;
        string id;

        while(!file.eof() && reading_init==true) {
            getline(file,line);
            id = line_get_first_word(line);
            if (id==rang_id) {
                int range = stoi(line_get_arg(line,2));
                if (range<0){
                    throw std::domain_error("IDB::read_file(): Invalid range specified. range cannot be negative!");
                }
                interaction_range = (unsigned int)range;
                rang_read = true;
            }
            if (id==flav_id) {
                bp_types = line_merge_tail(line,2);
                flav_read = true;
            }
            if (id==disc_id) {
                disc_len = stod(line_get_arg(line,2));
                disc_read = true;
            }
            if (id=="type") {
                break;
            }
        }
        if (!rang_read){
            throw std::domain_error("IDB::read_file(): Inconsistent Database File: Missing coupling range \n requires line: range = coupling_range");
        }
        if (!flav_read){
            throw std::domain_error("IDB::read_file(): Inconsistent Database File: Missing bp flavor \n requires line: flavors = string_of_bp_flavors");
        }
        if (!disc_read){
            throw std::domain_error("IDB::read_file(): Inconsistent Database File: Missing discretization length \n requires line: discretization = disc_val");
        }

        amount_units  = pow(bp_types.length(),2*(interaction_range+1));
        sets_per_unit = 1+2*interaction_range;

        units = new IDBUnit[amount_units];

        unsigned int unit_counter = 0;
        while(!file.eof()) {
            id = line_get_first_word(line);
            if (id=="type") {


                arma::cube Ms = arma::zeros(3,3,sets_per_unit);
                arma::colvec Theta0;
                string type = line_get_arg(line,1);

                for (unsigned int i=0;i<sets_per_unit;++i) {
                    getline(file,line);
                    id = line_get_first_word(line);
                    if (id != stiff_id) {
                        throw std::domain_error("IDB::read_file(): Inconsistent Database File: Insufficient matrices defined.");
                    }
                    Ms.slice(i) = line_get_mat(line,1);
                }
                getline(file,line);
                id = line_get_first_word(line);
                if (id != ground_id) {
                    throw std::domain_error("IDB::read_file(): Inconsistent Database File: Groundstate rotation vector definition missing.");
                }
                Theta0 = line_get_vec(line,1);

                if (IDBU_add_unit(type, Theta0, Ms,(unit_counter))) {
                    ++unit_counter;
                }
            }
            getline(file,line);
        }

        // TODO: offload the status prints to dedicated function
        cout << "Interaction range set to " << interaction_range << endl;
        cout << "Monomer types:   " << bp_types << endl;
        cout << "Discretization length set to "   << disc_len << endl;

        if (unit_counter!=amount_units) {
            cout << "Warning: IDB::read_file(): Missing interaction in Interaction Database" << endl;
            if (unit_counter==1) {
                cout << unit_counter << " Interaction Unit initialized (" << amount_units << " Units expected.)" << endl;
            }
            else {
                cout << unit_counter << " Interaction Units initialized (" << amount_units << " Units expected.)" << endl;
            }
        }
        else {
            if (amount_units==1) {
                cout << unit_counter << " Interaction Unit initialized" << endl;
            }
            else {
                cout << unit_counter << " Interaction Units initialized" << endl;
            }
        }

        if (sets_per_unit==1) {
            cout << sets_per_unit << " interaction per Unit " << endl;
        }
        else {
            cout << sets_per_unit << " interactions per Unit " << endl;
        }
        return true;

    }
    else {
        throw std::domain_error("IDB::read_file(): Interaction File not found!");
        return false;
    }
}

string IDB::line_get_first_word(const string& line) {
    vector<string> entries;
    split_string(line, entries);
    if (entries.size() == 0) {
        return "";
    }
    return entries[0];
}

string IDB::line_get_arg(const string& line,unsigned int num) {
    vector<string> entries;
    split_string(line, entries);
    if (entries.size() < (num+1)) {
        return "";
    }
    return entries[num];
}

string IDB::line_merge_tail(const string& line,unsigned int num) {
    vector<string> entries;
    split_string(line, entries);
    if (entries.size() < (num+1)) {
        return "";
    }
    string str = "";
    for (unsigned i=num;i<entries.size();++i) {
        str = str + entries[i];
    }
    return str;
}

arma::mat IDB::line_get_mat(const string& line,unsigned start) {
    vector<string> entries;
    split_string(line, entries);
    if (entries.size() != (start+9)) {
        cout << "Invalid matrix line:" << endl;
        cout << line << endl;
        throw std::domain_error("IDB::line_get_mat: Inconsistent Database File: invalid amount of entries in matrix line");
    }
    unsigned row;
    unsigned col;
    unsigned matid;
    arma::mat matrix = arma::zeros(3,3);
    for (unsigned i=start;i<(start+9);++i) {
        matid = i-start;
        col = matid%3;
        row = (matid-col)/3;
        matrix(row,col) = stod(entries[i]);
    }
    return matrix;
}

arma::colvec IDB::line_get_vec(const string& line,unsigned int start) {
    vector<string> entries;
    split_string(line, entries);
    if (entries.size() != (start+3)) {
        cout << "Invalid vector line:" << endl;
        cout << line << endl;
        throw std::domain_error("IDB::line_get_mat: Inconsistent Database File: invalid amount of entries in vector line");
    }
    unsigned int row;
    arma::colvec v = arma::zeros(3);
    for (unsigned i=start;i<(start+3);++i) {
        row = i-start;
        v(row) = stod(entries[i]);
    }
    return v;
}

void IDB::split_string(const string &str, vector<string> &output)
{
    string::size_type start = 0;
    string::size_type last;
    string::size_type spc_last;
    string::size_type tab_last;
    spc_last = str.find_first_of(" ");
    tab_last = str.find_first_of("\t");
    last = min(spc_last,tab_last);
    while (last != string::npos)
    {
        if (last > start)
        {
            output.push_back(str.substr(start, last - start));
        }
        start = ++last;
        spc_last = str.find_first_of(" ",last);
        tab_last = str.find_first_of("\t",last);
        last = min(spc_last,tab_last);
    }
    string last_str = str.substr(start);
    if (last_str.length() != 0) {
        if (last_str.at(0) != ' ' && last_str.at(0) != '\t') {
            output.push_back(str.substr(start));
        }
    }
}


