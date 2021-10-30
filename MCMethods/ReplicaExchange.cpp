#include "ReplicaExchange.h"


ReplicaExchange::ReplicaExchange(std::string inputfn, const std::vector<std::string> & argv_original) {
    std::vector<std::string> argv = argv_original;

    input = new InputRead(inputfn);

    // Check if all necessary arguments are contained in the input file
    if (!input->contains_multi("replica_exchange")) {
        throw std::invalid_argument("Error: ReplicaExchange: replica_exchange not specified in input file!");
    }
    MultiLineInstance * ML = input->get_multiline("replica_exchange")->get_instance(0);
    if (!ML->contains_singleline("temps")) {
        throw std::invalid_argument("Error: ReplicaExchange: argument 'temps' missing in multiline replica_exchange!");
    }
    if (!ML->contains_singleline("sweeps")) {
        throw std::invalid_argument("Error: ReplicaExchange: argument 'sweeps' missing in multiline replica_exchange!");
    }
    if (!ML->contains_singleline("sweepsteps")) {
        throw std::invalid_argument("Error: ReplicaExchange: argument 'sweepsteps' missing in multiline replica_exchange!");
    }

    // Initialize simulation paramters
    temps       = ML->get_single_vec("temps");
    sweeps      = ML->get_single_val<long long>("sweeps");
    sweepsteps  = ML->get_single_val<long long>("sweepsteps");
    if (ML->contains_singleline("output_all")) {
        output_all  = ML->get_single_val<bool>("output_all");
    }
    if (ML->contains_singleline("equi")) {
        equi  = ML->get_single_val<long long>("equi");
    }

    // Initialize statedump file. This file contains information about the current energy and temperature of
    // each replica
    std::string main_dump_dir;
    main_dump_dir = InputChoice_get_single<std::string> ("dump_dir",input,argv,DEFAULT_DUMP_DIR);
    main_dump_dir = InputChoice_get_single<std::string> ("dir",input,argv,main_dump_dir);
    statedump_filename = main_dump_dir + ".replica_exchange";
    std::ofstream ofstr;
    ofstr.open(statedump_filename, std::ofstream::out | std::ofstream::trunc);
    ofstr.close();

    // Initialize the seed
    seedstr = InputChoice_get_single<std::string>  ("seed",input,argv,seedstr);
    seedseq = seedstr2seedseq(seedstr);
    seedstr = seedseq2seedstr(seedseq);


    // Set seed for Replica Exchange random number generator
    std::seed_seq seed(seedseq.begin(), seedseq.end());
    gen.seed(seed);

    // Initialize main simulation (lowest energy)
    std::vector<std::string> main_argv;
    main_argv = add2argv(argv     ,"-seed",seedstr);
    main_argv = add2argv(main_argv,"-T",std::to_string(temps[0]));
    PolyMC * main_polyMC = new PolyMC(main_argv);
    std::cout << "Main Argv: " << std::endl;
    for (int j=0;j<main_argv.size();j++) {
        std::cout << main_argv[j] << " ";
    }

    polymcs .push_back(main_polyMC);
    energies.push_back(main_polyMC->get_chain()->extract_true_energy());
    betas   .push_back(main_polyMC->get_chain()->get_beta());
    sim_id  .push_back(0);

    // Initalize replica simulations
    for (int i=1;i<temps.size();i++) {
        std::cout << temps[i] << std::endl;
        std::cout << "----------------" << std::endl;

        std::vector<std::string> rep_argv;

//        nargv = add2argv(nargv,"-T",std::to_string(temps[i]));
        rep_argv = add2argv(argv,"-T",std::to_string(temps[i]));

        std::vector<long long int> nseedseq;
        for (int s=0;s<seedseq.size();s++) {
            nseedseq.push_back(seedseq[s]+i*1111);
        }
        std::string nseedstr = seedseq2seedstr(nseedseq);
        rep_argv = add2argv(rep_argv,"-seed",nseedstr);

        std::string dump_dir = main_dump_dir + "_" + std::to_string(temps[i]);
        rep_argv = add2argv(rep_argv,"-dump_dir",dump_dir);

        std::cout << "Replica Arg: " << std::endl;
        for (int j=0;j<rep_argv.size();j++) {
            std::cout << rep_argv[j] << " ";
        }


        PolyMC * new_polyMC = new PolyMC(rep_argv,output_all);
        polymcs .push_back(new_polyMC);
        energies.push_back(new_polyMC->get_chain()->extract_true_energy());
        betas   .push_back(new_polyMC->get_chain()->get_beta());
        sim_id  .push_back(i);
    }

    // Set range for exchange selection generator
    decltype(select_exchange.param()) new_range(0, polymcs.size()-1);
    select_exchange.param(new_range);

//    for (int i=0;i<1000;i++) {
//        replica_swap(0);
////        std::cout << select_exchange(gen) << std::endl;
//    }
//    std::cout << polymcs.size() <<  std::endl;
//    std::cout << std::endl;
//    std::exit(0);



    // equilibrate
    if (equi>0) {
        for (int i=0;i<polymcs.size();i++) {
//            std::string runname = "Equilibration T = " + std::to_string(temps[i]);
            polymcs[i]->run(equi,"Equilibration T = " + std::to_string(temps[i]),false,false);
        }
    }

    run();
}

ReplicaExchange::~ReplicaExchange() {
    delete input;
    for (unsigned i=0;i<polymcs.size();i++) {
        delete polymcs[i];
    }
}

bool ReplicaExchange::run() {
    for (int i=0;i<polymcs.size();i++) {
        polymcs[i]->init_external_run();
    }

    int print_every = 1000/sweepsteps;
    if (print_every == 0) {
        print_every = 1;
    }

    init_print_state();
    for (long long sweep=0;sweep<sweeps;sweep++) {
        for (int i=0;i<polymcs.size();i++) {
//            polymcs[i]->run(sweepsteps,"ReplicaExchange sweep " + std::to_string(sweep+1) +  " (T = " + std::to_string(temps[i])+")",true,false);
            polymcs[i]->external_run(sweepsteps,"ReplicaExchange sweep " + std::to_string(sweep+1) +  " (T = " + std::to_string(temps[i])+")",true,false,false,false,false);
        }
        replica_swap(sweep);

        if (sweep%print_every==0) {
            get_all_energies();
            dump_stats(sweep);
        }
        if (sweep%100==0) {
            print_state(sweep,sweeps);
        }

    }
    return true;
}

bool ReplicaExchange::replica_swap(long long int sweep) {
    double boltzmann;
    double randval;
    double tmpd;
    int    tmpi;
    int    id1,id2;
    int    updown;

    id1 = select_exchange(gen);
    if (id1 == 0) {
        id2 = 1;
    }
    else if (id1 == polymcs.size()-1) {
        id2 = id1 - 1;
    }
    else {
        updown = randbinary(gen);
        if (updown == 0) {
            updown = -1;
        }
        id2 = id1 + updown;
    }
    boltzmann = metropolis(polymcs[id1],polymcs[id2]);
    randval   = uniformdist(gen);
//    std::cout << boltzmann << std::endl;
    if (randval<boltzmann) {
        std::cout << "Swap " << id1 << " - " << id2 << std::endl;
        polymcs[id1]->exchange_configs(polymcs[id2]);

        tmpd          = energies[id1];
        energies[id1] = energies[id2];
        energies[id2] = tmpd;

        tmpi        = sim_id[id1];
        sim_id[id1] = sim_id[id2];
        sim_id[id2] = tmpi;

//        std::cout << energies[id1] << " " << polymcs[id1]->get_chain()->extract_true_energy() << std::endl;
//        std::cout << energies[id2] << " " << polymcs[id2]->get_chain()->extract_true_energy() << std::endl;
        return true;
    }
//    std::cout << id1 << " " << id2 << std::endl;
    return true;
}



int ReplicaExchange::replica_swaps(long long int sweep) {
    double boltzmann;
    double randval;
    double tmpd;
    int    tmpi;
    int    id1,id2;
    get_all_energies();

    int num_swaps = 0;

    std::vector<int> order = gen_rand_order(polymcs.size()-1);
    for (int i=0;i<polymcs.size()-1;i++) {
        id1 = order[i];
        id2 = id1+1;
        boltzmann = metropolis(polymcs[id1],polymcs[id2]);
        randval   = uniformdist(gen);
//        std::cout << boltzmann << std::endl;
        if (randval<boltzmann) {
//            std::cout << "Swap " << id1 << " - " << id2 << std::endl;
            polymcs[id1]->exchange_configs(polymcs[id2]);

            tmpd          = energies[id1];
            energies[id1] = energies[id2];
            energies[id2] = tmpd;

            tmpi        = sim_id[id1];
            sim_id[id1] = sim_id[id2];
            sim_id[id2] = tmpi;

            num_swaps++;
//            std::cout << energies[id1] << " " << polymcs[id1]->get_chain()->extract_true_energy() << std::endl;
//            std::cout << energies[id2] << " " << polymcs[id2]->get_chain()->extract_true_energy() << std::endl;
        }
    }
    return num_swaps;
}

std::vector<int> ReplicaExchange::gen_rand_order(int num) {
    std::vector<int> order;
    std::vector<int> nums;
    for (int i=0;i<num;i++) {
        nums.push_back(i);
    }
    for (int i=0;i<num;i++) {
        decltype(genorder.param()) new_range(0, num-i-1);
        genorder.param(new_range);
        int id = genorder(gen);
        order.push_back(nums[id]);
        nums.erase(nums.begin()+id);
    }
    return order;
}

double ReplicaExchange::metropolis(PolyMC* RunA, PolyMC* RunB) {
    double EA,EB,betaA,betaB;
    double dE;

    EA = RunA->get_chain()->extract_true_energy();
    EB = RunB->get_chain()->extract_true_energy();
    betaA = RunA->get_chain()->get_beta();
    betaB = RunB->get_chain()->get_beta();

    dE = (EA-EB)*(betaA-betaB);
    return std::exp(dE);
}

void ReplicaExchange::get_all_energies() {
    for (int i=0;i<energies.size();i++) {
        energies[i] = polymcs[i]->get_chain()->extract_true_energy();
    }
}

void ReplicaExchange::dump_stats(long long sweep) {
    std::ofstream ofstr;
    ofstr.open(statedump_filename, std::ofstream::out | std::ofstream::app);
    ofstr << sweep;
    for (int i=0;i<energies.size();i++) {
         ofstr << " " << energies[i];
    }
    for (int i=0;i<sim_id.size();i++) {
         ofstr << " " << sim_id[i];
    }
    ofstr << "\n";
    ofstr.close();
}

void ReplicaExchange::init_print_state() {
    timer_start = std::chrono::high_resolution_clock::now();
}

void ReplicaExchange::print_state(long long sweep,long long tot_sweeps) {
    timer_finish   = std::chrono::high_resolution_clock::now();
    timer_elapsed  = timer_finish-timer_start;
    std::cout << "#############################" << std::endl;
    std::cout << "sweep " << sweep << " (" << tot_sweeps << ")" << std::endl;
    std::cout << "  elapsed time:             " << time_remaining((timer_elapsed.count())/sweep*(sweep)) << "\n";
    std::cout << "  estimated to finished in: " << time_remaining((timer_elapsed.count())/sweep*(tot_sweeps-sweep)) << "\n";
}


std::string ReplicaExchange::time_remaining(int seconds) {
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

