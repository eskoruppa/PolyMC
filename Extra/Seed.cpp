#include "Seed.h"

std::vector<long long> seedstr2seedseq(const std::string & seedstr, int num_entries) {
    std::vector<long long> seedseq;
    if (seedstr == "-1") {
        seedseq = gen_seedseq(num_entries);
    }
    else {
        size_t pos = 0;
        std::string num;
        std::string cp_seedstr = seedstr;
        while ((pos = cp_seedstr.find("-")) != std::string::npos) {
            num = cp_seedstr.substr(0, pos);
            seedseq.push_back(std::stoll(num));
            cp_seedstr.erase(0, pos + 1);
        }
        seedseq.push_back(std::stoll(cp_seedstr));
    }
    return seedseq;
}

std::string seedseq2seedstr(const std::vector<long long> & seedseq) {
    std::string seedstr;
    std::ostringstream ss;
    ss << seedseq[0];
    seedstr = ss.str();
    for (int i=1;i<SEED_SEQ_LEN;i++) {
        std::ostringstream ss1;
        ss1 << seedseq[i];
        seedstr += "-"+ss1.str();
    }
    return seedstr;
}

std::vector<long long int> gen_seedseq(int num_entries) {
    std::vector<long long> seedseq;
    std::random_device  rd{};
    for (int i=0;i<num_entries;i++) {
        seedseq.push_back(rd());
    }
    return seedseq;
}






