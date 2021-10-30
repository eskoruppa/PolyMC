#include "Plasmid_old.h"

void plasmid() {

//cout.precision(17);

    string inter_fn = "IDB/Plasmid6bp";


    double T = 50;

    double sigma  = -0.00;
    int num_bp    = 100;

    int print  = 10000;
    int cal_Lk = 1;
    bool Lk_first = false;

    bool EV_active = false;
    bool kill_if_link_voilated = false;


    int turns = num_bp/10*(1.+sigma);
//    turns = 11;
    cout << "turns = " << turns << endl;

    Chain chain(inter_fn);
    chain.gen_circular(num_bp, sigma, "a");


    chain.set_T0_subtract(false);
    chain.set_T(T);


////// EXCLUDED VOLUMES /////

    double bp_per_EV  = 11;
    bp_per_EV = 8;
    double EV_rad     = 0.34*(bp_per_EV+0.1);
    int num_EV_per_PC = 11;


    ExVol EV(&chain,EV_rad,num_EV_per_PC);

    int rep_size_min = 1;
    int rep_size_max = 1;
    int rep_dist_min = 20;
    int rep_dist_max = 100;
    bool translate_twist = true;
    MCS_RepTherm rep(&chain,rep_size_min,rep_size_max,rep_dist_min,rep_dist_max,translate_twist);


    MCS_CSrot   rot1(&chain,2,40);
    MCS_CSrot   rot2(&chain,larger(2,num_bp/32),larger(2,num_bp/16));
    MCS_CSrot   rot3(&chain,larger(2,num_bp/16),larger(2,num_bp/8));
    MCS_CSrot   rot4(&chain,larger(2,num_bp/8),larger(2,num_bp/4));
    MCS_CSrot   rot5(&chain,larger(2,num_bp/4),larger(2,num_bp/2));

    MCS_Pivot   piv(&chain);
    MCS_CStrans trans(&chain,2,100,{0,0,1});
    MCS_Twist   twi(&chain);



    if (EV_active) {
        rep.set_excluded_volume(&EV);
        piv.set_excluded_volume(&EV);
        trans.set_excluded_volume(&EV);

        rot1.set_excluded_volume(&EV);
        rot2.set_excluded_volume(&EV);
        rot3.set_excluded_volume(&EV);
        rot4.set_excluded_volume(&EV);
        rot5.set_excluded_volume(&EV);
    }
    chain.fix_link();



    long long int steps = 1e10;
    long long int equi  = 1e4;


    std::ostringstream strs;
    strs << T;

    int    XYZ_dump_every      = 1000;
    string XYZ_dumpxyzfilename = "dump/plasmid"+strs.str()+".xyz";
    string XYZ_dumpxyzcenter   = "COM";
    string XYZ_rep             = "simple"; //"EV"

    int    triads_dump_every      = 1000;
    string triads_dumpxyzfilename = "dump/kink_circ.triads";
    string triads_dumpxyzcenter   = "COM";

    int    stiff_dump_every  = 50;
    string stiff_filename    = "dump/stiff";
    bool   stiff_all         = true;

    int    lb_dump_every  = 25;
    string lb_filename    = "dump/lb";
    int    lb_m_max       = 10;

    int    FE_dump_every  = 25;
    string FE_filename    = "dump/FE";

    int    EN_dump_every  = 10;
    string EN_filename    = "dump/Energy"+strs.str();

    int    ETS_dump_every  = 50;
    string ETS_filename    = "dump/Ceff";
    string ETS_writhe      = "fuller";

    int    LINK_dump_every  = 100;
    string LINK_filename    = "dump/Link";
    string LINK_options     = "quick";

    int    DPEP_dump_every         = 10000;
//        string DPEP_filename           = "dump/Branches/plasmid_endpoints_4.6k_#14";
    string DPEP_filename           = "dump/Branches/zzplasmid_endpoints_4.0k_#1";
    double DPEP_density            = 0.1;
    int    DPEP_sample_step_dist   = 1;
    int    DPEP_sample_points      = 1;
    bool   DPEP_dump_writhe        = true;
    bool   DPEP_dump_xyz           = true;


    int    DM_dump_every  = 50000;
    string DM_filename    = "dump/branch_test/distmap";
    double DM_density     = 0.1; //0.2;


    ///// DUMPS

    Dump_avgThetas          DTheta(&chain, 1e5, "dump/Thetas",false);
    Dump_Stiff              Dstiff(&chain, stiff_dump_every, stiff_filename,stiff_all);
    Dump_xyz                Dxyz(&chain, XYZ_dump_every, XYZ_dumpxyzfilename ,false, XYZ_dumpxyzcenter,XYZ_rep);
    Dump_xyz                Dtriads(&chain, triads_dump_every, triads_dumpxyzfilename ,false, triads_dumpxyzcenter,"pdb");
    Dump_PersLen            Dlb(&chain, lb_dump_every, lb_filename,true,lb_m_max);
    Dump_ForceExtension     DFE(&chain, FE_dump_every, FE_filename,0.25*num_bp,0.75*num_bp,true);
    Dump_Energy             DEN(&chain, EN_dump_every, EN_filename,false);
    Dump_EffTorsStiff       DETS(&chain, ETS_dump_every, ETS_filename, {0,0,1},ETS_writhe ,0.25*num_bp,0.75*num_bp,true);

    Dump_Linkingnumber      DLINK(&chain, LINK_dump_every, LINK_filename,LINK_options,false);
    Dump_DistMap            DDM(&chain, DM_dump_every, DM_filename, DM_density);
    Dump_PlasmidEndpoints   DPEP(&chain, DPEP_dump_every, DPEP_filename, DPEP_density,DPEP_sample_step_dist,DPEP_sample_points,DPEP_dump_writhe,DPEP_dump_xyz);


    double LK;
    double TW;

    double LK1_F = 0;
    double LK2_F = 0;
    double LK1 = 0;
    double LK2 = 0;
    int LK_count = 0;
    int LK_N = 1000;

    arma::colvec dir = {0,0,1};


    bool acpt;

    vector<BPStep*> BPS = *chain.get_BPS();

    auto timer_start = std::chrono::high_resolution_clock::now();
    auto timer_finish = std::chrono::high_resolution_clock::now();
    chrono::duration<double> timer_elapsed;

    for (long long int step=0;step<equi;step++) {
        if (step%10000==0) {
            cout << "Equilibration-step " << step << endl;
        }
//            acpt = piv.MC();
//            acpt = trans.MC();
        acpt = rot1.MC();
        acpt = rot2.MC();
        acpt = rot3.MC();
        acpt = rot4.MC();
        acpt = rot5.MC();
        Dxyz.dump();
    }
    cout << "Equilibrated" << endl;


    double Wr,Tw,Lk,fixed_LK;
    Wr = chain.cal_langowski_writhe_1a(0.2);
    Tw = chain.cal_twist(0,num_bp);
    fixed_LK = Wr+Tw;

    double dLk=0;


    for (long long int step=0;step<=steps;step++) {

//            Wr = chain.cal_langowski_writhe_1a(0.05);

        if (step%print==0) {

            if (step==50000) {
                chain.set_sigma(0);
                fixed_LK=0;
            }

            cout << "################" << endl;
            cout << "step " << step << endl;
            cout << "rep  acceptance rate = " << rep.acceptance_rate()   << endl;
            cout << "piv  acceptance rate = " << piv.acceptance_rate()   << endl;
            cout << "twi  acceptance rate = " << twi.acceptance_rate()   << endl;
            cout << "rot1 acceptance rate = " << rot1.acceptance_rate()   << endl;
            cout << "rot2 acceptance rate = " << rot2.acceptance_rate()   << endl;
            cout << "rot3 acceptance rate = " << rot3.acceptance_rate()   << endl;
            cout << "rot4 acceptance rate = " << rot4.acceptance_rate()   << endl;
            cout << "rot5 acceptance rate = " << rot5.acceptance_rate()   << endl;
            cout << "traj acceptance rate = " << trans.acceptance_rate() << endl;
            cout << "ExVol rejection rate = " << EV.rejection_rate() << endl;
//                cout << "Lang1a Writhe: " << chain.cal_langowski_writhe_1a() << endl;
//                cout << "Fuller Writhe: " << chain.cal_fuller_writhe(0,num_bp,dir) << endl;
//                cout << "Fuller Twist : " << chain.cal_twist(0,num_bp) << endl;

            timer_finish = std::chrono::high_resolution_clock::now();
            timer_elapsed       = timer_finish - timer_start;
            cout << "elapsed time:   " << timer_elapsed.count() << " s\n";

            timer_start = std::chrono::high_resolution_clock::now();

            if (!chain.config_consistent()) {
                cout << "chain positions inconsistent!" << endl;
                std::exit(0);
            }

            if (step%(cal_Lk*print)==0) {



                Wr = chain.cal_langowski_writhe_1a(1);
                Tw = chain.cal_twist();
                Lk = Wr+Tw;
                dLk = std::abs(Lk-fixed_LK);

                cout << "Wr = " << Wr << endl;
                cout << "Tw = " << Tw << endl;
                cout << "Lk = " << Lk << endl;

                cout << "Error langowski = " << Lk - fixed_LK << endl;

//                double gWr = chain.gauss_writhe(1);
//                cout << "Lk = " << Tw+gWr << endl;
//                cout << "Error gauss     = " << Tw+gWr - fixed_LK << endl;

                if (dLk > 1) {
                    if (kill_if_link_voilated) {
                        Dxyz.final_dump();
                        std::exit(0);
                    }
                }
            }
        }


//            for (int rep=0;rep<30;rep++) {
//                acpt = twi.MC();
//            }




//        if (step%1==0) {
//            acpt = rep.MC();
//            if (acpt) {
//                if (!chain.config_consistent()) {
//                    cout << "chain positions inconsistent!" << endl;
//                    cout << "reptation!" << endl;
//                    Dxyz.final_dump();
//                    std::exit(0);
//                }
//            }
//        }

        for (int r=0;r<5;r++) {
            acpt = rot1.MC();
        }

        acpt = rot2.MC();
        acpt = rot3.MC();
        acpt = rot4.MC();
        acpt = rot5.MC();





//            DTheta.dump();
        Dxyz.dump();
//        DEN.dump();
        Dtriads.dump();
//            Dstiff.dump();
//            Dlb.dump();
//            DFE.dump();
//            DETS.dump();
//            DDM.dump();
//        DPEP.dump();
//            DLINK.dump();


    }
//        DTheta.final_dump();
//        Dstiff.final_dump();
//        Dlb.final_dump();
//        DFE.final_dump();
//        DETS.final_dump();
//        DDM.final_dump();
//        DPEP.final_dump();

    cout << "piv acceptance rate = " << piv.acceptance_rate() << endl;
    cout << "rot acceptance rate = " << rot1.acceptance_rate() << endl;
    cout << "tra acceptance rate = " << trans.acceptance_rate() << endl;

}
