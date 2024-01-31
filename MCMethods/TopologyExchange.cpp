#include "TopologyExchange.h"


TopologyExchange::TopologyExchange(std::string inputfn, const std::vector<std::string> & argv_original) 
: swap_counter(0)
{
    std::vector<std::string> argv = argv_original;

    input = new InputRead(inputfn);

    // Check if all necessary arguments are contained in the input file
    if (!input->contains_multi("topology_exchange")) {
        throw std::invalid_argument("Error: TopologyExchange: topology_exchange not specified in input file!");
    }
    MultiLineInstance * ML = input->get_multiline("topology_exchange")->get_instance(0);
    if (!ML->contains_singleline("Lks")) {
        throw std::invalid_argument("Error: TopologyExchange: argument 'Lks' missing in multiline topology_exchange!");
    }
    if (!ML->contains_singleline("sweeps")) {
        throw std::invalid_argument("Error: TopologyExchange: argument 'sweeps' missing in multiline topology_exchange!");
    }
    if (!ML->contains_singleline("sweepsteps")) {
        throw std::invalid_argument("Error: TopologyExchange: argument 'sweepsteps' missing in multiline topology_exchange!");
    }


    // Initialize simulation paramters
    lks         = ML->get_single_vec("Lks");
    sweeps      = ML->get_single_val<long long>("sweeps");
    sweepsteps  = ML->get_single_val<long long>("sweepsteps");
    if (ML->contains_singleline("output_all")) {
        output_all  = ML->get_single_val<bool>("output_all");
    }
    if (ML->contains_singleline("equi")) {
        equi  = ML->get_single_val<long long>("equi");
    }
    else {
        equi = 0;
    }
    if (ML->contains_singleline("swaps")) {
        swap_attempts = ML->get_single_val<int>("swaps");
    }
    else {
        swap_attempts = lks.size() / 2;
    }

    if (!ML->contains_singleline("num_threads")) 
    {
        num_threads = 1;
    }
    else {
        num_threads  = ML->get_single_val<int>("num_threads");
    }
    if (!ML->contains_singleline("stats_every")) 
    {
        stats_every = 1;
    }
    else {
        stats_every  = ML->get_single_val<long long int>("stats_every");
    }


    use_openmp = num_threads > 1;

    // Initialize Dump dirs
    base_dump_dir = InputChoice_get_single<std::string> ("dump_dir",input,argv,DEFAULT_DUMP_DIR);
    base_dump_dir = InputChoice_get_single<std::string> ("dir",input,argv,base_dump_dir);

    // Initialize the seed
    seedstr = InputChoice_get_single<std::string>  ("seed",input,argv,seedstr);
    seedseq = seedstr2seedseq(seedstr);
    seedstr = seedseq2seedstr(seedseq);

    // Set seed for Replica Exchange random number generator
    std::seed_seq seed(seedseq.begin(), seedseq.end());
    gen.seed(seed);

    // Dummy for dlk to simga conversion
    PolyMC * dummy = new PolyMC(argv,false);

    // Initialize simulations
    for (int i=0;i<lks.size();i++) {
        std::cout << "Initializing PolyMC for LK=" << lks[i] << std::endl;

        double sigma = dummy->get_chain()->dlk2sigma(lks[i]);

        // Set new argument vector
        std::vector<std::string> rep_argv;
        
        //sigma
        std::stringstream sigstream;
        sigstream << std::fixed << std::setprecision(16) << sigma;
        std::string sigmastr = sigstream.str();
        rep_argv = add2argv(argv,"-sigma",sigmastr);

        // seed
        rep_argv = add2argv(rep_argv,"-seed",seedstr);

        // individual filename
        std::stringstream lkstream;
        lkstream << std::fixed << std::setprecision(2) << lks[i];
        std::string lkstr = lkstream.str();
        std::replace( lkstr.begin(), lkstr.end(), '.', 'p');
        std::string dump_dir = base_dump_dir + "_lk" + lkstr;
        rep_argv = add2argv(rep_argv,"-dump_dir",dump_dir);

        std::cout << "Replica Args: " << std::endl;
        for (int j=0;j<rep_argv.size();j++) {
            std::cout << rep_argv[j] << " ";
        }

        PolyMC * new_polymc = new PolyMC(rep_argv,output_all);
        // new_polymc->get_chain()->set_Delta_Lk(lks[i]);
        // new_polymc->set_all_backups();
        // std::cout << "dlk = " << new_polymc->get_chain()->get_dLK() << std::endl;

        Replica * new_replica = new Replica(new_polymc, "replica_"+lkstr);
        new_replica->dprops["dlk"]    = lks[i];
        new_replica->dprops["sigma"]  = sigma;
        new_replica->iprops["id"]     = i;
        new_replica->iprops["confid"] = i;
        new_replica->dprops["energy"] = new_polymc->get_chain()->extract_full_energy();
        
        replicas.push_back(new_replica);
    }

    // delete dummy PolyMC
    delete dummy;

    // init stats dump
    init_stats();

    // Set range for exchange selection generator
    decltype(select_exchange.param()) new_range(0, replicas.size()-2);
    select_exchange.param(new_range);


    #ifdef TOPOLEXCHANGE_OPENMP
    // Set number of threads
    std::cout << "number of threads = " << num_threads << std::endl;
    omp_set_num_threads(num_threads);
    #endif

    // equilibrate
    if (equi>0) {
        #ifdef TOPOLEXCHANGE_OPENMP
        #pragma omp parallel for
        #endif
        for (int i=0;i<replicas.size();i++) {
            std::string runname = "Equilibration lk = " + std::to_string(lks[i]);
            std::cout << "Equilibration of " << runname << " (" << equi << " steps)" << std::endl;

            replicas[i]->add_steps(equi);
            replicas[i]->start_timer();
            replicas[i]->polymc->run(equi,"Equilibration lk = " + std::to_string(lks[i]),false,false);
            replicas[i]->stop_timer();
        }
    }
    run();
}

TopologyExchange::~TopologyExchange() {
    delete input;
    for (unsigned i=0;i<replicas.size();i++) {
        delete replicas[i]->polymc;
        delete replicas[i];
    }
}

bool TopologyExchange::run() {

    for (int i=0;i<replicas.size();i++) {
        replicas[i]->polymc->init_external_run();
    }
    init_print_state();
    dump_stats(0);

    for (long long sweep=1;sweep<=sweeps;sweep++) { 

        double mean_rate = mean_steprate();
        #ifdef TOPOLEXCHANGE_OPENMP
        #pragma omp parallel for
        #endif
        for (int i=0;i<replicas.size();i++) {
            long long indiv_sweepsteps = replicas[i]->steprate() / mean_rate * sweepsteps;
            // std::cout << "lk " << replicas[i]->dprops["dlk"] << " - <rate> = " << replicas[i]->steprate() << " (" << indiv_sweepsteps << " steps)" << std::endl;
            /*
                TODO: round invid_sweepsteps to next largest multiple of largest dumpn or to value passed in topology_exchange block
            */

            replicas[i]->add_steps(sweepsteps);
            replicas[i]->start_timer();
            replicas[i]->polymc->external_run(sweepsteps,"TopologyExchange sweep " + std::to_string(sweep+1) +  " (T = " + std::to_string(lks[i])+")",true,false,false,false,false);
            replicas[i]->stop_timer();
        }

        // std::cout << "mean rate = " << mean_rate << " steps/s" << std::endl;
        // for (int i=0;i<replicas.size();i++) {
        //     long long indiv_sweepsteps = replicas[i]->steprate() / mean_rate * sweepsteps;
        //     std::cout << "lk " << replicas[i]->dprops["dlk"] << " - <rate> = " << replicas[i]->steprate() << " (" << indiv_sweepsteps << " steps)" << std::endl;
        // }

        replica_swaps(sweep,swap_attempts);

        if (sweep%stats_every==0) {
            dump_stats(sweep);
            print_state(sweep,sweeps);
        }
    }
    return true;
}

double TopologyExchange::mean_steprate() {
    double crate = 0;
    for (int i=0;i<replicas.size();i++) {
        crate += replicas[i]->steprate();
    }
    return crate / replicas.size();
}

bool TopologyExchange::replica_swaps(long long int sweep, int num_swaps) {

    std::vector<int> swapids;
    // generate swap ids
    for (int i=0;i<num_swaps;i++) {
        swapids.push_back(select_exchange(gen));
    }
    
    // std::cout << "##############" << std::endl;
    // for (int j=0; j<swapids.size(); j++) {
    //     std::cout << swapids.at(j) << ' ';
    // }
    // std::cout << "\n";

    #ifdef TOPOLEXCHANGE_OPENMPSWAPS
    int curr = 0;
    std::vector<std::vector<int>> partials;    
    while (curr < num_swaps) {
        std::vector<int> partial;
        int first = curr;
        int last  = curr;
        partial.push_back(swapids[curr]);
        while (curr < num_swaps) {
            curr += 1;
            int nextid = swapids[curr];
            bool sufficientspacing = true;
            for (unsigned i=0;i<partial.size();i++) {
                if (std::abs(partial[i]-nextid) < 2) {
                    sufficientspacing = false;
                    break;
                }
            }
            if (sufficientspacing) {
                partial.push_back(nextid);
            }
            else {
                partials.push_back(partial);
                break;
            }
        }
    }

    // for (int i=0;i<partials.size();i++) {
    //     for (int j=0; j<partials.at(i).size(); j++) {
    //         std::cout << partials.at(i).at(j) << ' ';
    //     }
    //     std::cout << "\n";
    // }

    for (int i=0;i<partials.size();i++) {
        if (partials.at(i).size() > 1) {
            #ifdef TOPOLEXCHANGE_OPENMP
            #pragma omp parallel for
            #endif
            for (int j=0; j<partials.at(i).size(); j++) {
                int id1 = partials.at(i).at(j);
                int id2 = id1+1;
                if (replica_swap(sweep,id1,id2)) {
                    swap_counter += 1;
                }
            }
        }
        else {
            for (int j=0; j<partials.at(i).size(); j++) {
                int id1 = partials.at(i).at(j);
                int id2 = id1+1;
                if (replica_swap(sweep,id1,id2)) {
                    swap_counter += 1;
                }
            }
        }
    }
    #else
    for (int i=0;i<swapids.size();i++) {
        int id1 = swapids.at(i);
        int id2 = id1+1;
        if (replica_swap(sweep,id1,id2)) {
            swap_counter += 1;
        }
    }
    #endif

    return true;
}


bool TopologyExchange::replica_swap(long long int sweep, int id1, int id2) {
    double E1,E1p;
    double E2,E2p;
    double dE1,dE2;
    double boltzmann;
    double randval;

    // // TESTING
    // std::cout << "##############################" << std::endl;
    // int nbp = replicas[id1]->polymc->get_chain()->get_num_bp();
    // std::cout << nbp-1 << std::endl;
    // std::cout << id1 << " " << id2 << std::endl;
    // std::cout << replicas[id1]->polymc->get_chain()->get_triads()->slice(nbp-1) << std::endl;
    // std::cout << replicas[id2]->polymc->get_chain()->get_triads()->slice(nbp-1) << std::endl;
    // ///////////////////////////////

    E1 = replicas[id1]->polymc->get_chain()->extract_full_energy();
    E2 = replicas[id2]->polymc->get_chain()->extract_full_energy();

    replicas[id1]->polymc->get_chain()->set_Delta_Lk(replicas[id2]->dprops["dlk"]);
    replicas[id2]->polymc->get_chain()->set_Delta_Lk(replicas[id1]->dprops["dlk"]);
    replicas[id1]->polymc->exchange_configs(replicas[id2]->polymc);

    // std::cout << replicas[id1]->dprops["dlk"] << " " <<  replicas[id1]->polymc->cal_dlk() << std::endl;
    // std::cout << replicas[id2]->dprops["dlk"] << " " <<  replicas[id2]->polymc->cal_dlk() << std::endl;

    E1p = replicas[id1]->polymc->get_chain()->extract_full_energy();
    E2p = replicas[id2]->polymc->get_chain()->extract_full_energy();

    dE1 = E1p-E1;
    dE2 = E2p-E2;
    double beta = replicas[id1]->polymc->get_chain()->get_beta();
    beta = 1;
    boltzmann = std::exp(-beta*(dE1+dE2));
    randval   = uniformdist(gen);

    if (randval<boltzmann) {
        // std::cout << "swapped " << id1 << " and " << id2 << std::endl;
        // std::cout << "dE1 = " << dE1 << std::endl;
        // std::cout << "dE2 = " << dE2 << std::endl;
        // std::cout << "p = " << boltzmann << std::endl;

        replicas[id1]->polymc->set_all_backups();
        replicas[id2]->polymc->set_all_backups();

        int tmp_confid = replicas[id1]->iprops["confid"];
        replicas[id1]->iprops["confid"] = replicas[id2]->iprops["confid"];
        replicas[id2]->iprops["confid"] = tmp_confid;
        return true;
    }

    replicas[id1]->polymc->get_chain()->set_Delta_Lk(replicas[id2]->dprops["dlk"]);
    replicas[id2]->polymc->get_chain()->set_Delta_Lk(replicas[id1]->dprops["dlk"]);
    replicas[id1]->polymc->exchange_configs(replicas[id2]->polymc);
    replicas[id1]->polymc->get_chain()->extract_full_energy();
    replicas[id2]->polymc->get_chain()->extract_full_energy();

    // std::cout << replicas[id1]->polymc->get_chain()->get_triads()->slice(nbp-1) << std::endl;
    // std::cout << replicas[id2]->polymc->get_chain()->get_triads()->slice(nbp-1) << std::endl;
    // std::exit(0);

    return false;
}

void TopologyExchange::init_stats() {
    replicastate_filename = base_dump_dir + ".replica_exchange";
    std::ofstream ofstr;
    ofstr.open(replicastate_filename, std::ofstream::out | std::ofstream::trunc);
    ofstr << sweeps;
    for (int i=0;i<replicas.size();i++) {
         ofstr << " " << replicas[i]->dprops["dlk"];
    }
    ofstr << "\n";
    ofstr.close();
}

void TopologyExchange::dump_stats(long long sweep) {
    std::ofstream ofstr;
    ofstr.open(replicastate_filename, std::ofstream::out | std::ofstream::app);
    ofstr << sweep;
    for (int i=0;i<replicas.size();i++) {
         ofstr << " " << replicas[i]->iprops["confid"];
    }
    ofstr << "\n";
    ofstr.close();
}

void TopologyExchange::init_print_state() {
    timer_start = std::chrono::high_resolution_clock::now();
}

void TopologyExchange::print_state(long long sweep,long long tot_sweeps) {
    timer_finish   = std::chrono::high_resolution_clock::now();
    timer_elapsed  = timer_finish-timer_start;
    std::cout << "#############################" << std::endl;
    std::cout << "sweep " << sweep << " (" << tot_sweeps << ")" << std::endl;
    std::cout << "number of swaps:            " << swap_counter << std::endl;
    std::cout << "  elapsed time:             " << time_remaining((timer_elapsed.count())/sweep*(sweep)) << "\n";
    std::cout << "  estimated to finished in: " << time_remaining((timer_elapsed.count())/sweep*(tot_sweeps-sweep)) << "\n";
}


std::string TopologyExchange::time_remaining(int seconds) {
    int days    = seconds/86400;
    seconds     = seconds%86400;
    int hours   = seconds/3600;
    seconds     = seconds%3600;
    int minutes = seconds/60;
    seconds     = seconds%60;
    std::string str;
    if (days>0) {
        str += std::to_string(days) + "d ";
    }
    if (hours>0) {
        str += std::to_string(hours) + "h ";
    }
    if (minutes>0) {
        str += std::to_string(minutes) + "min ";
    }
    str += std::to_string(seconds) + "s ";
    return str;
}
