#include "Replica.h"


Replica::Replica(PolyMC * pmc, const std::string & name) 
: polymc(pmc), name(name), step_count(0)
{
    elapsed = 0;
    timer_running = false;
}

Replica::~Replica() {
}

void Replica::add_dprop(const std::string name, double prop) {
    dprops[name] = prop;
}

void Replica::add_iprop(const std::string name, int prop) {
    iprops[name] = prop;
}

void Replica::add_steps(int steps) {
    step_count += steps;
}

void Replica::reset_steps() {
    step_count = 0;
}

void Replica::start_timer() {
    if (timer_running) {
        throw std::runtime_error("Attempting to start timer while timer was still running.");
    }
    timer_running = true;
    timer_start = std::chrono::high_resolution_clock::now();
}

void Replica::stop_timer() {
    if (!timer_running) {
        throw std::runtime_error("Attempting to stop timer while timer was not running.");
    }
    timer_finish = std::chrono::high_resolution_clock::now();
    timer_running = false;
    timer_elapsed = timer_finish - timer_start;
    elapsed += timer_elapsed.count();
}

void Replica::reset_elapsed_time() {
    elapsed = 0;
}

double Replica::steprate() {
    return step_count / elapsed;
}