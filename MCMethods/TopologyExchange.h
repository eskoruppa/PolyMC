#ifndef __TOPOLOGY_EXCHANGE_INCLUDED__
#define __TOPOLOGY_EXCHANGE_INCLUDED__

#include <omp.h>
#include <unistd.h>
#include <chrono>

#include "../Input/Argparse.h"
#include "../Input/InputRead.h"
#include "../Input/ManipulateArgv.h"
#include "../PolyMC.h"
#include "Replica.h"

#define TOPOLEXCHANGE_OPENMP
#define TOPOLEXCHANGE_OPENMPSWAPS


class TopologyExchange;
class TopologyExchange {

protected:

    long long sweeps     = 0;
    long long sweepsteps = 0;
    int swap_attempts    = 1;
    long long equi       = 0;
    bool output_all      = false;
    bool use_openmp      = false;
    int  num_threads = 1;
    int swap_counter;

    std::string                 seedstr = "-1";
    std::vector<long long int>  seedseq;

    std::vector<Replica*> replicas;


    std::vector<PolyMC*> polymcs;
    std::vector<double>  lks;
    // std::vector<double>  sigmas;
    // std::vector<double>  energies;
    // std::vector<double>  betas;
    // std::vector<int>     sim_id;

    std::string base_dump_dir;
    std::string replicastate_filename;

    InputRead * input;

    std::mt19937                            gen{0};
    std::uniform_real_distribution<double>  uniformdist{0.0,1.0};
    std::uniform_int_distribution<>         genorder{0, 1};

    std::uniform_int_distribution<>         select_exchange{0, 1};
    std::uniform_int_distribution<>         randbinary{0, 1};

    std::chrono::high_resolution_clock::time_point timer_start;
    std::chrono::high_resolution_clock::time_point timer_finish;
    std::chrono::high_resolution_clock::time_point timer_ref;
    std::chrono::duration<double> timer_elapsed;

    long long int stats_every; 


public:
    TopologyExchange(std::string input_file, const std::vector<std::string> & argv);
    ~TopologyExchange();

    bool run();

    double mean_steprate();
    void init_stats();
    void dump_stats(long long sweep);

protected:
    bool replica_swaps(long long int sweep, int num_swaps);
    bool replica_swap(long long int sweep, int id1, int id2);
    

    void        init_print_state();
    void        print_state(long long sweep,long long tot_sweeps);
    std::string time_remaining(int seconds);

};








#endif
