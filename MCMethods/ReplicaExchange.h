#ifndef __REPLICA_EXCHANGE_INCLUDED__
#define __REPLICA_EXCHANGE_INCLUDED__
#include "../Input/Argparse.h"
#include "../Input/InputRead.h"
#include "../Input/ManipulateArgv.h"
#include "../PolyMC.h"



class ReplicaExchange;

class ReplicaExchange {

protected:

    long long sweeps     = 0;
    long long sweepsteps = 0;
    long long equi       = 0;
    bool output_all      = false;

    std::string                 seedstr = "-1";
    std::vector<long long int>  seedseq;

    std::vector<double>  temps;
    std::vector<PolyMC*> polymcs;
    std::vector<double>  energies;
    std::vector<double>  betas;
    std::vector<int>     sim_id;

    std::string statedump_filename;

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


public:
    ReplicaExchange(std::string input_file, const std::vector<std::string> & argv);
    ~ReplicaExchange();

    bool run();

    void dump_stats(long long sweep);

protected:
    int    replica_swaps(long long int sweep);
    bool   replica_swap(long long int sweep);
    std::vector<int> gen_rand_order(int num);
    double metropolis(PolyMC* RunA, PolyMC* RunB);
    void   get_all_energies();

    void        init_print_state();
    void        print_state(long long sweep,long long tot_sweeps);
    std::string time_remaining(int seconds);

};






#endif
