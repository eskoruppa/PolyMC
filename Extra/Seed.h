#ifndef __SEED_INCLUDED__
#define __SEED_INCLUDED__

//// include input functions
#include "../Input/Argparse.h"
#include "../Input/InputRead.h"
#include "../Input/InputChoice.h"
#include "../Input/GenInputFile.h"

#include <chrono>
#include <armadillo>
#include <cmath>
#include <vector>
#include <random>

#define SEED_SEQ_LEN 4

std::vector<long long int> seedstr2seedseq(const std::string & seedstr, int num_entries=SEED_SEQ_LEN);
std::string seedseq2seedstr(const std::vector<long long int> & seedseq);
std::vector<long long int> gen_seedseq(int num_entries=SEED_SEQ_LEN);

#endif

