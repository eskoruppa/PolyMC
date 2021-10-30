#include "WriteLoadRestart.h"



std::vector<RestartData> loadrestart(std::string filename){

    std::vector<RestartData> datavec;

    std::ifstream file;
    std::string   line;
    std::vector<std::string> linelist;
    std::vector<std::string> lines;

    file.open(filename);
    if (file.is_open())
    {
        while(!file.eof()) {
            std::getline(file,line);
            lines.push_back(line);
        }
    }
    else {
        std::cout << "Error: Restart file not found!" << std::endl;
        std::exit(1);
    }

    unsigned lid = 0;
    while (lid<lines.size()-1) {
        if (lines[lid].at(0) == '#') {
            lid++;
            get_linelist(lines[lid],linelist);
            if (linelist[0] != "snapshot:") {
                continue;
            }
            long long int snapshot = std::stoll(linelist[1]);

            lid++;
            get_linelist(lines[lid],linelist);
            if (linelist[0] != "type:") {
                continue;
            }
            std::string type = linelist[1];

            lid++;
            get_linelist(lines[lid],linelist);
            if (linelist[0] != "num_bp:") {
                continue;
            }
            long long int num_bp = std::stoll(linelist[1]);

            lid++;
            get_linelist(lines[lid],linelist);
            if (linelist[0] != "sequence:") {
                continue;
            }
            std::string seq = linelist[1];

            lid++;
            get_linelist(lines[lid],linelist);
            if (linelist[0] != "dLK:") {
                continue;
            }
            double dLK = std::stod(linelist[1]);

            if (lines.size()-1-lid < num_bp) {
                break;
            }

            arma::mat pos = arma::zeros(3,num_bp);
            for (unsigned i=0;i<num_bp;i++) {
                lid++;
                get_linelist(lines[lid],linelist);
                pos.col(i)(0) = std::stod(linelist[0]);
                pos.col(i)(1) = std::stod(linelist[1]);
                pos.col(i)(2) = std::stod(linelist[2]);
            }
            arma::cube triads = arma::zeros(3,3,num_bp);
            for (unsigned i=0;i<num_bp;i++) {
                lid++;
                get_linelist(lines[lid],linelist);
                triads.slice(i).col(0)(0) = std::stod(linelist[0]);
                triads.slice(i).col(0)(1) = std::stod(linelist[1]);
                triads.slice(i).col(0)(2) = std::stod(linelist[2]);
                triads.slice(i).col(1)(0) = std::stod(linelist[3]);
                triads.slice(i).col(1)(1) = std::stod(linelist[4]);
                triads.slice(i).col(1)(2) = std::stod(linelist[5]);
                triads.slice(i).col(2)(0) = std::stod(linelist[6]);
                triads.slice(i).col(2)(1) = std::stod(linelist[7]);
                triads.slice(i).col(2)(2) = std::stod(linelist[8]);
            }

            RestartData restart_data;
            restart_data.triads   = triads;
            restart_data.pos      = pos;
            restart_data.num_bp   = num_bp;
            restart_data.type     = type;
            restart_data.sequence = seq;
            restart_data.snapshot = snapshot;
            restart_data.dLK      = dLK;
            datavec.push_back(restart_data);

        }
        lid++;
    }
    return datavec;
}

