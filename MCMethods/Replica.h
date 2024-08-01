#ifndef __REPLICA_INCLUDED__
#define __REPLICA_INCLUDED__

#include <chrono>

#include <unordered_map>
#include "../Input/Argparse.h"
#include "../Input/InputRead.h"
#include "../Input/ManipulateArgv.h"
#include "../PolyMC.h"


class Replica;

class Replica {

public:

    PolyMC * polymc;
    std::string name;
    double  energy;

    std::string dump_dir;
    bool output_all = true;
    std::unordered_map<std::string, double> dprops;
    std::unordered_map<std::string, int>    iprops;

    long long step_count;

    std::chrono::high_resolution_clock::time_point timer_start;
    std::chrono::high_resolution_clock::time_point timer_finish;
    std::chrono::high_resolution_clock::time_point timer_ref;
    std::chrono::duration<double> timer_elapsed;
    double elapsed;
    bool timer_running;

public:

    Replica(PolyMC * pmc, const std::string & name);
    ~Replica();

    void add_dprop(const std::string name, double prop);
    void add_iprop(const std::string name, int prop);

    void add_steps(int steps);
    void reset_steps();

    void start_timer();
    void stop_timer();
    void reset_elapsed_time();

    double steprate();

};

#endif