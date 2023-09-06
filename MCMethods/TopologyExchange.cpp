#include "TopologyExchange.h"


TopologyExchange::TopologyExchange(std::string inputfn, const std::vector<std::string> & argv_original) {
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
    use_openmp = num_threads > 1;

    // Initialize Dump dirs
    base_dump_dir = InputChoice_get_single<std::string> ("dump_dir",input,argv,DEFAULT_DUMP_DIR);
    base_dump_dir = InputChoice_get_single<std::string> ("dir",input,argv,base_dump_dir);
    replicastate_filename = base_dump_dir + ".replica_exchange";
    std::ofstream ofstr;
    ofstr.open(replicastate_filename, std::ofstream::out | std::ofstream::trunc);
    ofstr.close();

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
        std::stringstream stream;

        //sigma
        stream << std::fixed << std::setprecision(16) << sigma;
        std::string sigmastr = stream.str();
        rep_argv = add2argv(argv,"-sigma",sigmastr);

        // seed
        rep_argv = add2argv(rep_argv,"-seed",seedstr);

        // individual filename
        stream << std::fixed << std::setprecision(2) << lks[i];
        std::string lkstr = stream.str();
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

    for (long long sweep=0;sweep<sweeps;sweep++) { 
        std::cout << "sweep " << sweep << std::endl;
        
        double mean_rate = mean_steprate();

        std::cout << "mean rate = " << mean_rate << " steps/s" << std::endl;

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

        for (int i=0;i<replicas.size();i++) {
            long long indiv_sweepsteps = replicas[i]->steprate() / mean_rate * sweepsteps;
            std::cout << "lk " << replicas[i]->dprops["dlk"] << " - <rate> = " << replicas[i]->steprate() << " (" << indiv_sweepsteps << " steps)" << std::endl;
        }

        for (int s=0;s<=swap_attempts;s++) {
            replica_swap(sweep);
        }

        dump_stats(sweep);

    }

//     int print_every = 1000/sweepsteps;
//     if (print_every == 0) {
//         print_every = 1;
//     }

//     init_print_state();
//     for (long long sweep=0;sweep<sweeps;sweep++) {
//         for (int i=0;i<polymcs.size();i++) {
// //            polymcs[i]->run(sweepsteps,"TopologyExchange sweep " + std::to_string(sweep+1) +  " (T = " + std::to_string(lks[i])+")",true,false);
//             polymcs[i]->external_run(sweepsteps,"TopologyExchange sweep " + std::to_string(sweep+1) +  " (T = " + std::to_string(lks[i])+")",true,false,false,false,false);
//         }
//         replica_swap(sweep);

//         if (sweep%print_every==0) {
//             get_all_energies();
//             dump_stats(sweep);
//         }
//         if (sweep%100==0) {
//             print_state(sweep,sweeps);
//         }

//     }
    return true;
}

double TopologyExchange::mean_steprate() {
    double crate = 0;
    for (int i=0;i<replicas.size();i++) {
        crate += replicas[i]->steprate();
    }
    return crate / replicas.size();
}



bool TopologyExchange::replica_swap(long long int sweep) {
    int id1,id2;
    double E1,E1p;
    double E2,E2p;
    double dE1,dE2;
    double boltzmann;
    double randval;

    id1 = select_exchange(gen);
    id2 = id1 + 1;
    // std::cout << "id1 = " << id1 << std::endl;

    E1 = replicas[id1]->polymc->get_chain()->extract_full_energy();
    E2 = replicas[id2]->polymc->get_chain()->extract_full_energy();

    // std::cout << "swap attempt" << std::endl;
    // std::cout << replicas[id1]->dprops["dlk"] << " " <<  replicas[id1]->polymc->cal_dlk() << std::endl;
    // std::cout << replicas[id2]->dprops["dlk"] << " " <<  replicas[id2]->polymc->cal_dlk() << std::endl;

    replicas[id1]->polymc->get_chain()->set_Delta_Lk(replicas[id2]->dprops["dlk"]);
    replicas[id2]->polymc->get_chain()->set_Delta_Lk(replicas[id1]->dprops["dlk"]);
    replicas[id1]->polymc->exchange_configs(replicas[id2]->polymc);

    // std::cout << replicas[id1]->dprops["dlk"] << " " <<  replicas[id1]->polymc->cal_dlk() << std::endl;
    // std::cout << replicas[id2]->dprops["dlk"] << " " <<  replicas[id2]->polymc->cal_dlk() << std::endl;

    E1p = replicas[id1]->polymc->get_chain()->extract_full_energy();
    E2p = replicas[id2]->polymc->get_chain()->extract_full_energy();

    dE1 = E1p-E1;
    dE2 = E2p-E2;
    boltzmann = std::exp(-dE1-dE2);
    randval   = uniformdist(gen);

    if (randval<boltzmann) {
        std::cout << "swapped " << id1 << " and " << id2 << std::endl;
        std::cout << "dE1 = " << dE1 << std::endl;
        std::cout << "dE2 = " << dE2 << std::endl;
        std::cout << "p = " << boltzmann << std::endl;

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
    return false;
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

// void TopologyExchange::init_print_state() {
//     timer_start = std::chrono::high_resolution_clock::now();
// }

// void TopologyExchange::print_state(long long sweep,long long tot_sweeps) {
//     timer_finish   = std::chrono::high_resolution_clock::now();
//     timer_elapsed  = timer_finish-timer_start;
//     std::cout << "#############################" << std::endl;
//     std::cout << "sweep " << sweep << " (" << tot_sweeps << ")" << std::endl;
//     std::cout << "  elapsed time:             " << time_remaining((timer_elapsed.count())/sweep*(sweep)) << "\n";
//     std::cout << "  estimated to finished in: " << time_remaining((timer_elapsed.count())/sweep*(tot_sweeps-sweep)) << "\n";
// }


// std::string TopologyExchange::time_remaining(int seconds) {
//     int days    = seconds/86400;
//     seconds     = seconds%86400;
//     int hours   = seconds/3600;
//     seconds     = seconds%3600;
//     int minutes = seconds/60;
//     seconds     = seconds%60;
//     std::string str;
//     if (days>0) {
//         str += std::to_string(days) + "d ";
//     }
//     if (hours>0) {
//         str += std::to_string(hours) + "h ";
//     }
//     if (minutes>0) {
//         str += std::to_string(minutes) + "min ";
//     }
//     str += std::to_string(seconds) + "s ";
//     return str;
// }
