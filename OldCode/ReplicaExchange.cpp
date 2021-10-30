#include "Plasmid.h"

void plasmid_replica_exchange() {

//cout.precision(17);

    string interaction_file = "IDB/Plasmid6bp";
    string type             = "closed";

    int T = 300;
    vector<double> Temps;
    Temps.push_back(T);

    double factor     = 1.08;
    int    n_replicas = 13;

    for (int i=0;i<n_replicas;i++) {
        T *= 1.08;
        Temps.push_back(T);
    }

//    Temps.push_back(300);
//    Temps.push_back(330);
//    Temps.push_back(365);
//    Temps.push_back(400);
//    Temps.push_back(440);
//    Temps.push_back(485);
//    Temps.push_back(535);
//    Temps.push_back(590);
//    Temps.push_back(650);

    long long int sweeps     = 2e6;
    long long int steps      = 50;
    long long int equi_steps = 1e4;

    double sigma    = -0.05;
    int num_bp      = 737;
    string sequence = "a";

    double EV_rad     = 3.5;
    int num_EV_per_PB = 11;

    int dump_every = 200;
    int dump_every_replica = 5000;

//    string dump_dir = "dump/Branches/LongRuns/3.5nm_7.37kbp_#3";
    string dump_dir = "dump/Branches/test/2600";
    dump_dir += "/Run";
    string dump_fn;

    std::ostringstream strs;

    vector<PolyMC*> Runs;

    cout << "#############################################" << endl;
    cout << "#############################################" << endl;
    cout << "Initializing Production Replica" << endl;
    PolyMC Run1;
    Run1.setup(type, interaction_file, num_bp, sequence, sigma, Temps[0]);
    cout << "PolyMC - setup done" << endl;
    Run1.set_ExVol(EV_rad,num_EV_per_PB, true);
    cout << "PolyMC - ExVol initialized" << endl;
    strs.clear();
    strs.str("");
    strs << Temps[0];
    dump_fn = dump_dir + "_" + strs.str() + "K_";
    Run1.set_dumps(dump_fn,dump_every);
    cout << "PolyMC - Dumps initialized" << endl;
    Runs.push_back(&Run1);

//    cout << "#############################################" << endl;
//    cout << "#############################################" << endl;
//    cout << "Initializing Replica 2" << endl;
//    PolyMC Run2;
//    Run2.setup(type, interaction_file, num_bp, sequence, sigma, Temps[1]);
//    cout << "PolyMC - setup done" << endl;
//    Run2.set_ExVol(EV_rad,num_EV_per_PB, true);
//    cout << "PolyMC - ExVol initialized" << endl;
//    strs.clear();
//    strs.str("");
//    strs << Temps[1];
//    dump_fn = dump_dir + "_" + strs.str() + "K_";
//    Run2.set_dumps(dump_fn);
//    cout << "PolyMC - Dumps initialized" << endl;
//    Runs.push_back(&Run2);

//
//    cout << "#############################################" << endl;
//    cout << "#############################################" << endl;
//    cout << "Initializing Replica 3" << endl;
//    PolyMC Run3;
//    Run3.setup(type, interaction_file, num_bp, sequence, sigma, Temps[2]);
//    cout << "PolyMC - setup done" << endl;
//    Run3.set_ExVol(EV_rad,num_EV_per_PB, true);
//    cout << "PolyMC - ExVol initialized" << endl;
//    strs.clear();
//    strs.str("");
//    strs << Temps[2];
//    dump_fn = dump_dir + "_" + strs.str() + "K_";
//    Run3.set_dumps(dump_fn);
//    cout << "PolyMC - Dumps initialized" << endl;
//    Runs.push_back(&Run3);
//
//    cout << "#############################################" << endl;
//    cout << "#############################################" << endl;
//    cout << "Initializing Replica 4" << endl;
//    PolyMC Run4;
//    Run4.setup(type, interaction_file, num_bp, sequence, sigma, Temps[3]);
//    cout << "PolyMC - setup done" << endl;
//    Run4.set_ExVol(EV_rad,num_EV_per_PB, true);
//    cout << "PolyMC - ExVol initialized" << endl;
//    strs.clear();
//    strs.str("");
//    strs << Temps[3];
//    dump_fn = dump_dir + "_" + strs.str() + "K_";
//    Run4.set_dumps(dump_fn);
//    cout << "PolyMC - Dumps initialized" << endl;
//    Runs.push_back(&Run4);

    for (int i=1;i<Temps.size();i++) {
        cout << "#############################################" << endl;
        cout << "#############################################" << endl;
        cout << "Initializing Replica " << i << endl;
        PolyMC * New_Run = new PolyMC;
        New_Run->setup(type, interaction_file, num_bp, sequence, sigma, Temps[i]);
        cout << "PolyMC - setup done" << endl;
        New_Run->set_ExVol(EV_rad,num_EV_per_PB, true);
        cout << "PolyMC - ExVol initialized" << endl;
        strs.clear();
        strs.str("");
        strs << Temps[i];
        dump_fn = dump_dir + "_" + strs.str() + "K_";
        New_Run->set_dumps(dump_fn,dump_every_replica);
        cout << "PolyMC - Dumps initialized" << endl;
        Runs.push_back(New_Run);
    }

    cout << endl << "#############################################" << endl;
    cout << "Equilibrating Simulations" << endl;
    for (int r=0;r<Runs.size();r++) {
        cout << "Equilibrating replica " << r+1 << " .. " << endl;
        Runs[r]->run(equi_steps,false, false);
        cout << "done" << endl;
    }
    cout << "#############################################" << endl;
    cout << "Commencing Production Run...\n" << endl;
    RE_run_SERIAL(Runs,Temps,steps,sweeps,dump_dir);
}


bool RE_run_SERIAL(vector<PolyMC*>& Runs,vector<double>& Temps, long long int steps, long long int sweeps, string& dump_dir) {

    vector<int> swap_counts;
    for (int i=0;i<Runs.size()-1;i++) {
        swap_counts.push_back(0);
    }

    vector<int> current_position;
    for (int i=0;i<Runs.size();i++) {
        current_position.push_back(i+1);
    }

    int cp_temp;
    string current_position_fn = dump_dir + "_current_positions";
    ofstream ofstr;
    ofstr.open(current_position_fn, ofstream::out | ofstream::trunc);
    ofstr.close();

    std::random_device                      rd{};
    std::mt19937                            gen{rd()};
    std::uniform_real_distribution<double>  uniformdist{0.0,1.0};
    std::seed_seq seed{rd(), rd(), rd(), rd(), rd(), rd(), rd(), rd()};
    gen.seed(seed);
    double boltzmann,randval;

    for (int sweep=0;sweep<sweeps;sweep++) {

        if (sweep%10==9) {
            cout << "###########################" << endl;
            cout << "sweep " << (sweep+1) << endl;
        }


        for (int r=0;r<Runs.size();r++) {
//            cout << "T = " << Temps[r] << endl;
            Runs[r]->run(steps);

            if (r>0) {
                boltzmann = ReplicaExchange_Metropolis(Runs[r-1],Runs[r]);
                randval = uniformdist(gen);
                if (randval<boltzmann) {
                    cout << "Swap " << r << " - " << r+1 << endl;
                    ReplicaExchange_Swap(Runs[r-1],Runs[r]);

                    cp_temp = current_position[r-1];
                    current_position[r-1] = current_position[r];
                    current_position[r] = cp_temp;

                    swap_counts[r-1]++;
                }
            }
        }

        ofstream ofstr;
        ofstr.open(current_position_fn, ofstream::app);
        for (int r=0;r<Runs.size();r++) {
            ofstr <<  current_position[r] << " ";
        }
        ofstr << "\n";
        ofstr.close();


    }
}





bool RE_run_OMP(vector<PolyMC*>& Runs,vector<double>& Temps, long long int steps, long long int sweeps, string& dump_dir) {

    vector<double> loadtime;
    vector<double> steptime;
    vector<int>    individual_steps;

    for (int i=0;i<Runs.size();i++) {
        loadtime.push_back(0);
        steptime.push_back(0);
        individual_steps.push_back(steps);
    }
    double time_per_step,min_time_per_step,min_time,max_time;

    vector<int> swap_counts;
    for (int i=0;i<Runs.size()-1;i++) {
        swap_counts.push_back(0);
    }

    omp_set_num_threads(Runs.size());
    for (int sweep=0;sweep<sweeps;sweep++) {

        #pragma omp parallel
        {
            #pragma omp for
            for (int proc=0;proc<Runs.size();proc++) {
                loadtime[proc] = Runs[proc]->run(individual_steps[proc]);
            }
        }

        cout << "###########################" << endl;
        cout << "sweep " << (sweep+1) << endl;
        cout << "Load Times:" << endl;
        min_time_per_step = loadtime[0];
        min_time = min_time_per_step;
        max_time = 0;
        for (int i=0;i<Runs.size();i++) {
            time_per_step = loadtime[i]/individual_steps[i];
            steptime[i]   = time_per_step;
            cout << "Run " << i+1 << ":  " << loadtime[i] << "s (" << individual_steps[i] << " steps - " << steptime[i] << "s per step)" << endl;
            if (time_per_step < min_time_per_step) {
                min_time_per_step = time_per_step;
            }

            if (loadtime[i] < min_time) {
                min_time = loadtime[i];
            }
            if (loadtime[i] > max_time) {
                max_time = loadtime[i];
            }
        }
        cout << "Load Wait      " << (max_time-loadtime[0])/loadtime[0]*100 << "%" << endl;
        cout << "Load Imbalance " << (max_time-min_time)/min_time*100 << "%" << endl;
//        for (int i=1;i<Runs.size();i++) {
//            if (steptime[i] > steptime[0]) {
//                individual_steps[i] = steps*steptime[0]/steptime[i];
//            }
//            else {
//                individual_steps[i] = steps;
//            }
//        }
//
//        ReplicaExchange(Runs,swap_counts);
    }
}


bool ReplicaExchange(vector<PolyMC*>& Runs, vector<int>& swap_counts) {
    int N = Runs.size();

    std::random_device                      rd{};
    std::mt19937                            gen{rd()};
    std::uniform_real_distribution<double>  uniformdist{0.0,1.0};
    std::seed_seq seed{rd(), rd(), rd(), rd(), rd(), rd(), rd(), rd()};
    gen.seed(seed);

    double boltzmann,randval;
    for (int i=0;i<N-1;i++) {
        boltzmann = ReplicaExchange_Metropolis(Runs[i],Runs[i+1]);
        randval = uniformdist(gen);
        if (randval<boltzmann) {
            cout << "Swap " << i << " - " << i+1 << endl;
            ReplicaExchange_Swap(Runs[i],Runs[i+1]);
            swap_counts[i]++;
        }
    }
    cout << "Swap counts: " << endl;
    for (int i=0;i<Runs.size()-1;i++) {
        cout << i << "-" << i+1 << " " << swap_counts[i] << endl;
    }
}

double ReplicaExchange_Metropolis(PolyMC* RunA, PolyMC* RunB) {
    double betaEA, betaEB, EA,EB,TA,TB;
    double dE;

    betaEA = RunA->get_chain()->extract_energy();
    betaEB = RunB->get_chain()->extract_energy();
    TA = RunA->get_chain()->get_T();
    TB = RunB->get_chain()->get_T();
    EA = betaEA*TA;
    EB = betaEB*TB;

    dE = (EA-EB)*(1/TA-1/TB);
    return std::exp(dE);
}


void ReplicaExchange_Swap(PolyMC* RunA, PolyMC* RunB) {
    arma::mat   temp_bp_pos = *RunA->get_chain()->get_bp_pos();
    arma::cube  temp_triads = *RunA->get_chain()->get_triads();

    RunA->get_chain()->set_config(RunB->get_chain()->get_bp_pos(),RunB->get_chain()->get_triads(),RunA->get_chain()->topology_closed());
    RunA->get_EV()->set_backup_conf(RunB->get_chain()->get_bp_pos(),RunB->get_chain()->get_triads());

    RunB->get_chain()->set_config(&temp_bp_pos,&temp_triads,RunB->get_chain()->topology_closed());
    RunB->get_EV()->set_backup_conf(&temp_bp_pos,&temp_triads);
}














