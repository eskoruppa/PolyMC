#include "InitDump.h"

std::vector<Dump*> init_dump_cmdargs(const std::vector<std::string> & argv, GenInputFile * geninfile, Chain * _chain, std::string mode, double EV_rad, std::string inputfn) {

    /*
    TODO: Add Dump Info

    bool print_info = false;
    if (mode == "dumpinfo" || mode == "dump_info") {
        print_info = true;
    }
    */


    int num_bp  = _chain->get_num_bp();
    std::string dump_dir = "";

    bool inputfile_given=false;
    InputRead input(inputfn);
    if (input.filefound()){
        dump_dir  = input.get_single_val<std::string>("dump_dir",dump_dir);
        inputfile_given=true;
    }
    dump_dir    = parse_arg(dump_dir, "-dir"     ,argv);
    dump_dir    = parse_arg(dump_dir, "-dump_dir",argv);

    // create path if it doesn't exist
    std::string path = std::filesystem::path(dump_dir).parent_path().u8string();
    // std::cout << path << std::endl;
    if (!fs::is_directory(path)) {
        fs::create_directories(path);
    }

    bool append_dumps = false;
    append_dumps = parse_flag(append_dumps,"-app", argv);
    append_dumps = InputChoice_get_single<bool>        ("app"     ,&input,argv,append_dumps);
    append_dumps = InputChoice_get_single<bool>        ("append"  ,&input,argv,append_dumps);

    std::vector<Dump*>   Dumps;

    /*
        Dump Meltingbubble Statistics
    */
    int    BUB_dump_every    = 0;
    std::string BUB_filename = dump_dir+".bub";
    BUB_dump_every  = InputChoice_get_single<int>        ("BUBn"  ,&input,argv,BUB_dump_every);
    BUB_dump_every  = InputChoice_get_single<int>        ("bubn"  ,&input,argv,BUB_dump_every);
    BUB_dump_every  = InputChoice_get_single<int>        ("Bubn"  ,&input,argv,BUB_dump_every);
    BUB_filename    = InputChoice_get_single<std::string>("BUBfn" ,&input,argv,BUB_filename);
    BUB_filename    = InputChoice_get_single<std::string>("bubfn" ,&input,argv,BUB_filename);
    BUB_filename    = InputChoice_get_single<std::string>("Bubfn" ,&input,argv,BUB_filename);

    if (BUB_dump_every>0) {
        Dump_BubbleStats* Dbub = new Dump_BubbleStats(_chain, BUB_dump_every, BUB_filename ,append_dumps);
        Dumps.push_back(Dbub);

        geninfile->add_entry(GENINFILE_DUMPS,"BUBn",BUB_dump_every);
        geninfile->add_entry(GENINFILE_DUMPS,"BUBfn",BUB_filename);
        geninfile->add_newline(GENINFILE_DUMPS);
    }

    /*
        Dump Configuration
    */
    int    CONF_dump_every    = 0;
    std::string CONF_filename = dump_dir+".conf";
    CONF_dump_every  = InputChoice_get_single<int>        ("CONFn"  ,&input,argv,CONF_dump_every);
    CONF_filename    = InputChoice_get_single<std::string>("CONFfn" ,&input,argv,CONF_filename);

    if (CONF_dump_every>0) {
        Dump_Conf* Dconf = new Dump_Conf(_chain, CONF_dump_every, CONF_filename ,append_dumps);
        Dumps.push_back(Dconf);

        geninfile->add_entry(GENINFILE_DUMPS,"CONFn",CONF_dump_every);
//        geninfile->add_entry(GENINFILE_DUMPS,"CONFfn",CONF_filename);
        geninfile->add_newline(GENINFILE_DUMPS);
    }


    /*
        Dump Distance Map
    */

    int    DM_dump_every    = 0;
    double DM_density       = 1;
    std::string DM_filename      = dump_dir+".dm";
    DM_dump_every  = InputChoice_get_single<int>        ("DMn"   ,&input,argv,DM_dump_every);
    DM_density     = InputChoice_get_single<double>     ("DMdens",&input,argv,DM_density);
    DM_filename    = InputChoice_get_single<std::string>("DMfn"  ,&input,argv,DM_filename);
    if (DM_dump_every>0) {
        Dump_DistMap* DDM = new Dump_DistMap(_chain, DM_dump_every, DM_filename, DM_density, append_dumps);
        Dumps.push_back(DDM);

        geninfile->add_entry(GENINFILE_DUMPS,"DMn",DM_dump_every);
        geninfile->add_entry(GENINFILE_DUMPS,"DMdens",DM_density);
//        geninfile->add_entry(GENINFILE_DUMPS,"DMfn",DM_filename);
        geninfile->add_newline(GENINFILE_DUMPS);
    }

    /*
        Dump Linking Number from End-Bead Tracing
    */
    int EBT_dump_every = 0;
    std::string EBT_filename = dump_dir+".endlink";
    EBT_dump_every  = InputChoice_get_single<int>          ("EBn" ,&input,argv,EBT_dump_every);
    EBT_filename    = InputChoice_get_single<std::string>  ("EBfn",&input,argv,EBT_filename);
    if (EBT_dump_every>0) {
        Dump_Endbead* EBT = new Dump_Endbead(_chain, EBT_dump_every, EBT_filename, append_dumps);
        Dumps.push_back(EBT);

        geninfile->add_entry(GENINFILE_DUMPS,"EBn",EBT_dump_every);
//        geninfile->add_entry(GENINFILE_DUMPS,"EBfn",EBT_filename);
        geninfile->add_newline(GENINFILE_DUMPS);
    }

    /*
        Dump First and Second Momment of Linking Number from End-Bead Tracing
    */
    int EBTSTAT_dump_every = 0;
    std::string EBTSTAT_filename = dump_dir+".endlinkstats";
    EBTSTAT_dump_every  = InputChoice_get_single<int>         ("EBSTATn" ,&input,argv,EBTSTAT_dump_every);
    EBTSTAT_filename    = InputChoice_get_single<std::string> ("EBSTATfn",&input,argv,EBTSTAT_filename);
    if (EBTSTAT_dump_every>0) {
        Dump_Endbead_Stats* EBTSTAT = new Dump_Endbead_Stats(_chain, EBTSTAT_dump_every, EBTSTAT_filename, append_dumps);
        Dumps.push_back(EBTSTAT);

        geninfile->add_entry(GENINFILE_DUMPS,"EBSTATn" ,EBTSTAT_dump_every);
//        geninfile->add_entry(GENINFILE_DUMPS,"EBSTATfn",EBTSTAT_filename);
        geninfile->add_newline(GENINFILE_DUMPS);
    }

    /*
        Dump Ceff from End-Bead Tracing
    */
    int EBCEFF_dump_every = 0;
    std::string EBCEFF_filename = dump_dir+".ebceff";
    EBCEFF_dump_every  = InputChoice_get_single<int>         ("EBCEFFn" ,&input,argv,EBCEFF_dump_every);
    EBCEFF_filename    = InputChoice_get_single<std::string> ("EBCEFFfn",&input,argv,EBCEFF_filename);
    if (EBCEFF_dump_every>0) {
        Dump_Endbead_Ceff * EBCEFF = new Dump_Endbead_Ceff(_chain, EBCEFF_dump_every, EBCEFF_filename, true);
        Dumps.push_back(EBCEFF);

        geninfile->add_entry(GENINFILE_DUMPS,"EBCEFFn" ,EBCEFF_dump_every);
//        geninfile->add_entry(GENINFILE_DUMPS,"EBCEFFfn",EBCEFF_filename);
        geninfile->add_newline(GENINFILE_DUMPS);
    }

    /*
        Dump End-to-end distance
    */

    int    E2E_dump_every   = 0;
    std::string E2E_filename     = dump_dir+".e2e";
    E2E_dump_every  = InputChoice_get_single<int>         ("E2En" ,&input,argv,E2E_dump_every);
    E2E_filename    = InputChoice_get_single<std::string> ("E2Efn",&input,argv,E2E_filename);

    if (E2E_dump_every>0) {
        Dump_EndToEndDistance* E2E = new Dump_EndToEndDistance(_chain, E2E_dump_every, E2E_filename, append_dumps);
        Dumps.push_back(E2E);

        geninfile->add_entry(GENINFILE_DUMPS,"E2En" ,E2E_dump_every);
//        geninfile->add_entry(GENINFILE_DUMPS,"E2Efn",E2E_filename);
        geninfile->add_newline(GENINFILE_DUMPS);
    }

    /*
        Dump Energy
    */
    int    Energy_dump_every      = 0;
    std::string Energy_filename   = dump_dir+".en";
    Energy_dump_every  = InputChoice_get_single<int>         ("En" ,&input,argv,Energy_dump_every);
    Energy_filename    = InputChoice_get_single<std::string> ("Efn",&input,argv,Energy_filename);
    if (Energy_dump_every>0) {
        Dump_Energy* DEnergy = new Dump_Energy(_chain, Energy_dump_every, Energy_filename ,false,append_dumps);
        Dumps.push_back(DEnergy);

        geninfile->add_entry(GENINFILE_DUMPS,"En" ,Energy_dump_every);
//        geninfile->add_entry(GENINFILE_DUMPS,"Efn",Energy_filename);
        geninfile->add_newline(GENINFILE_DUMPS);
    }

    /*
        Dump Extension
    */
    int         Ext_dump_every = 0;
    std::string Ext_filename   = dump_dir+".zext";
    int         Ext_dump_from  = -1;
    int         Ext_dump_to    = -1;
    int         Ext_subdomsize = 0;

    Ext_dump_every  = InputChoice_get_single<int>         ("Extn" ,&input,argv,Ext_dump_every);
    Ext_dump_every  = InputChoice_get_single<int>         ("extn" ,&input,argv,Ext_dump_every);
    Ext_dump_every  = InputChoice_get_single<int>         ("EXTn" ,&input,argv,Ext_dump_every);
    Ext_filename    = InputChoice_get_single<std::string> ("Extfn",&input,argv,Ext_filename);
    Ext_filename    = InputChoice_get_single<std::string> ("extfn",&input,argv,Ext_filename);
    Ext_filename    = InputChoice_get_single<std::string> ("EXTfn",&input,argv,Ext_filename);

    Ext_dump_from   = InputChoice_get_single<int>         ("Extfr" ,&input,argv,Ext_dump_from);
    Ext_dump_from   = InputChoice_get_single<int>         ("extfr" ,&input,argv,Ext_dump_from);
    Ext_dump_from   = InputChoice_get_single<int>         ("EXTfr" ,&input,argv,Ext_dump_from);
    Ext_dump_to     = InputChoice_get_single<int>         ("Extto" ,&input,argv,Ext_dump_to);
    Ext_dump_to     = InputChoice_get_single<int>         ("extto" ,&input,argv,Ext_dump_to);
    Ext_dump_to     = InputChoice_get_single<int>         ("EXTto" ,&input,argv,Ext_dump_to);

    Ext_subdomsize  = InputChoice_get_single<int>         ("Extsds" ,&input,argv,Ext_subdomsize);
    Ext_subdomsize  = InputChoice_get_single<int>         ("Extsub" ,&input,argv,Ext_subdomsize);
    Ext_subdomsize  = InputChoice_get_single<int>         ("extsds" ,&input,argv,Ext_subdomsize);
    Ext_subdomsize  = InputChoice_get_single<int>         ("extsub" ,&input,argv,Ext_subdomsize);
    Ext_subdomsize  = InputChoice_get_single<int>         ("EXTsds" ,&input,argv,Ext_subdomsize);
    Ext_subdomsize  = InputChoice_get_single<int>         ("EXTsub" ,&input,argv,Ext_subdomsize);

    if (Ext_dump_every>0) {
        geninfile->add_entry(GENINFILE_DUMPS,"Extn" ,Ext_dump_every);
//        geninfile->add_entry(GENINFILE_DUMPS,"Extfn",Ext_filename);
        geninfile->add_newline(GENINFILE_DUMPS);

        // Dump_Extension* DEXT = new Dump_Extension(_chain, Ext_dump_every, Ext_filename,append_dumps);
        // Dumps.push_back(DEXT);
    
        if (Ext_subdomsize == 0) {
            Dump_Extension* DEXT = new Dump_Extension(_chain, Ext_dump_every, Ext_filename, Ext_dump_from, Ext_dump_to, append_dumps);
            Dumps.push_back(DEXT);
        }
        else {
            Dump_Extension* DEXT = new Dump_Extension(_chain, Ext_dump_every, Ext_filename, Ext_dump_from, Ext_dump_to, Ext_subdomsize, append_dumps);
            Dumps.push_back(DEXT);
        }

    }

    /*
        Dump Force Extension
    */
    int    FE_dump_every = 0;
    double FE_frac       = 1;
    std::string FE_filename   = dump_dir+".fe";

    FE_dump_every  = InputChoice_get_single<int>         ("fen" ,&input,argv,FE_dump_every);
    FE_dump_every  = InputChoice_get_single<int>         ("FEn" ,&input,argv,FE_dump_every);
    FE_filename    = InputChoice_get_single<std::string> ("fefn",&input,argv,FE_filename);
    FE_filename    = InputChoice_get_single<std::string> ("FEfn",&input,argv,FE_filename);
    FE_frac        = InputChoice_get_single<double>      ("fefrac" ,&input,argv,FE_frac);
    FE_frac        = InputChoice_get_single<double>      ("FEfrac" ,&input,argv,FE_frac);

    if (FE_dump_every>0) {
        if (FE_frac > 1) FE_frac = 1;
        if (FE_frac < 0) FE_frac = 0;
        int start_id = (1-FE_frac)*0.5*num_bp;
        int end_id   = (0.5 + FE_frac*0.5) *num_bp;
        Dump_ForceExtension* DFEXT = new Dump_ForceExtension(_chain, FE_dump_every, FE_filename,start_id,end_id,append_dumps);
        Dumps.push_back(DFEXT);

        geninfile->add_entry(GENINFILE_DUMPS,"FEn"   ,FE_dump_every);
//        geninfile->add_entry(GENINFILE_DUMPS,"FEfn"  ,FE_filename);
        geninfile->add_entry(GENINFILE_DUMPS,"FEfrac",FE_frac);
        geninfile->add_newline(GENINFILE_DUMPS);
    }

    /*
        Dump Full Coupling Matrix
    */

    std::string FCM_filename = dump_dir;
    bool FCM_dump = false;
    FCM_dump      = parse_flag(false    , "-FCM",argv);
    FCM_dump      = parse_flag(FCM_dump , "-fcm",argv);
    if (!FCM_dump) {
        FCM_dump      = InputChoice_get_single<bool>        ("FCM"  ,&input,argv,FCM_dump);
        FCM_dump      = InputChoice_get_single<bool>        ("fcm"  ,&input,argv,FCM_dump);
    }
    FCM_filename  = InputChoice_get_single<std::string> ("fcmfn",&input,argv,FCM_filename);
    FCM_filename  = InputChoice_get_single<std::string> ("FCMfn",&input,argv,FCM_filename);

    if (FCM_dump) {
        Dump_Full_Coup_Matrix(_chain,FCM_filename);

//        Dump_Full_Coup_Matrix* DFCM = new Dump_Full_Coup_Matrix(_chain, FE_dump_every, FE_filename,start_id,end_id,append_dumps);
        geninfile->add_entry(GENINFILE_DUMPS,"FCM"   ,FCM_dump);
//        geninfile->add_entry(GENINFILE_DUMPS,"FCMfn" ,FCM_filename);
        geninfile->add_newline(GENINFILE_DUMPS);
    }

    /*
        Dump HatCurve
    */
    int    HC_dump_every    = 0;
    int    HC_dump_to_file  = 10000;
    std::string HC_filename = dump_dir+".hc";

    HC_dump_every   = InputChoice_get_single<int>         ("HCn" ,&input,argv,HC_dump_every);
    HC_dump_to_file = InputChoice_get_single<int>         ("HCn2f" ,&input,argv,HC_dump_to_file);
    HC_filename     = InputChoice_get_single<std::string> ("HCfn",&input,argv,HC_filename);

    if (HC_dump_every>0) {
        Dump_HatCurve* DHC  = new Dump_HatCurve(_chain, HC_dump_every, HC_filename,HC_dump_to_file,append_dumps);
        Dumps.push_back(DHC);

        geninfile->add_entry(GENINFILE_DUMPS,"HCn"   ,HC_dump_every);
        geninfile->add_entry(GENINFILE_DUMPS,"HCn2f" ,HC_dump_to_file);
//        geninfile->add_entry(GENINFILE_DUMPS,"HCfn" ,HC_filename);
        geninfile->add_newline(GENINFILE_DUMPS);
    }

    int    HCS_dump_every    = 0;
    std::string HCS_filename      = dump_dir+".hcs";
    HCS_dump_every = InputChoice_get_single<int>         ("HCSn" ,&input,argv,HCS_dump_every);
    HCS_filename   = InputChoice_get_single<std::string> ("HCSfn",&input,argv,HCS_filename);

    if (HCS_dump_every>0) {
        Dump_HatCurveStatistics* DHCS = new Dump_HatCurveStatistics(_chain, HCS_dump_every, HCS_filename,append_dumps);
        Dumps.push_back(DHCS);

        geninfile->add_entry(GENINFILE_DUMPS,"HCSn"  ,HCS_dump_every);
//        geninfile->add_entry(GENINFILE_DUMPS,"HCSfn" ,HCS_filename);
        geninfile->add_newline(GENINFILE_DUMPS);
    }

    /*
        Dump Helicity
    */
    int    HEL_dump_every    = 0;
    int    HEL_dump_to_file  = 0;
    std::string HEL_filename = dump_dir+".helicity";

    HEL_dump_every   = InputChoice_get_single<int>         ("HELn"   ,&input,argv,HEL_dump_every);
    HEL_dump_to_file = InputChoice_get_single<int>         ("HELn2f" ,&input,argv,HEL_dump_to_file);
    HEL_filename     = InputChoice_get_single<std::string> ("HELfn"  ,&input,argv,HEL_filename);

    if (HEL_dump_every>0) {
        Dump_Helicity* DHEL = new Dump_Helicity(_chain, HEL_dump_every,HEL_dump_to_file, HEL_filename, append_dumps);
        Dumps.push_back(DHEL);

        geninfile->add_entry(GENINFILE_DUMPS,"HELn"   ,HEL_dump_every);
        geninfile->add_entry(GENINFILE_DUMPS,"HELn2f" ,HEL_dump_to_file);
//        geninfile->add_entry(GENINFILE_DUMPS,"HELfn"  ,HEL_filename);
        geninfile->add_newline(GENINFILE_DUMPS);
    }


    /*
        Dump Kink statistics
    */
    int         kinkstat_dump_every = 0;
    std::string kinkstat_filename   = dump_dir+".kinkstat";

    kinkstat_dump_every = InputChoice_get_single<int>         ("kinkstatn"    ,&input,argv,kinkstat_dump_every);
    kinkstat_filename   = InputChoice_get_single<std::string> ("kinkstatfn"   ,&input,argv,kinkstat_filename);


    if (kinkstat_dump_every>0) {
        Dump_Kinkstats* Dkinkstat = new Dump_Kinkstats(_chain, kinkstat_dump_every, kinkstat_filename, append_dumps);
        Dumps.push_back(Dkinkstat);

        geninfile->add_entry(GENINFILE_DUMPS,"kinkstatn"    ,kinkstat_dump_every);
        geninfile->add_entry(GENINFILE_DUMPS,"kinkstatfn"   ,kinkstat_filename);
        geninfile->add_newline(GENINFILE_DUMPS);
    }

    /*
        Dump Linking Number
    */
    int    LK_dump_every      = 0;
    std::string LK_filename   = dump_dir+".lk";
    std::string LK_option     = "exact";
    double LK_frac            = 1;

    LK_dump_every = InputChoice_get_single<int>         ("LKn"    ,&input,argv,LK_dump_every);
    LK_filename   = InputChoice_get_single<std::string> ("LKfn"   ,&input,argv,LK_filename);
    LK_option     = InputChoice_get_single<std::string> ("LKoptn" ,&input,argv,LK_option);
    LK_frac       = InputChoice_get_single<double>      ("LKfrac" ,&input,argv,LK_frac);

    if (LK_dump_every>0) {
        if (LK_frac > 1) LK_frac = 1;
        if (LK_frac < 0) LK_frac = 0;
        int start_id  = (1-LK_frac)*0.5*num_bp;
        int end_id    = (0.5 + LK_frac*0.5) *num_bp;

        Dump_Linkingnumber* DLK = new Dump_Linkingnumber(_chain, LK_dump_every, LK_filename, start_id, end_id, LK_option, append_dumps);
        Dumps.push_back(DLK);

        geninfile->add_entry(GENINFILE_DUMPS,"LKn"    ,LK_dump_every);
//        geninfile->add_entry(GENINFILE_DUMPS,"LKfn"   ,LK_filename);
        geninfile->add_entry(GENINFILE_DUMPS,"LKoptn" ,LK_option);
        geninfile->add_newline(GENINFILE_DUMPS);
    }

    /*
        Dump MSD
    */
    int    MSD_dump_every       = 0;
    std::string MSD_filename    = dump_dir+".msd";
    int    MSD_dist             = 50;

    MSD_dump_every = InputChoice_get_single<int>         ("MSDn"    ,&input,argv,MSD_dump_every);
    MSD_filename   = InputChoice_get_single<std::string> ("MSDfn"   ,&input,argv,MSD_filename);
    MSD_dist       = InputChoice_get_single<int>         ("MSDdist" ,&input,argv,MSD_dist);

    if (MSD_dump_every>0) {
        Dump_MSD* DMSD = new Dump_MSD(_chain, MSD_dump_every, MSD_filename, MSD_dist, append_dumps);
        Dumps.push_back(DMSD);

        geninfile->add_entry(GENINFILE_DUMPS,"MSDn"    ,MSD_dump_every);
//        geninfile->add_entry(GENINFILE_DUMPS,"MSDfn"   ,MSD_filename);
        geninfile->add_entry(GENINFILE_DUMPS,"MSDdist" ,MSD_dist);
        geninfile->add_newline(GENINFILE_DUMPS);
    }

    /*
        Dump Persistence Length
    */
    int    PersLen_dump_every       = 0;
    std::string PersLen_filename    = dump_dir+".lb";
    int    PersLen_dist             = 50;
    bool   PersLen_cal_tors         = false;
    PersLen_dump_every = InputChoice_get_single<int>         ("LBn"    ,&input,argv,PersLen_dump_every);
    PersLen_filename   = InputChoice_get_single<std::string> ("LBfn"   ,&input,argv,PersLen_filename);
    PersLen_dist       = InputChoice_get_single<int>         ("LBdist" ,&input,argv,PersLen_dist);
    PersLen_cal_tors   = InputChoice_get_single<bool>        ("LBtors" ,&input,argv,PersLen_cal_tors);

    if (PersLen_dump_every>0) {
        Dump_PersLen* DLB = new Dump_PersLen(_chain, PersLen_dump_every, PersLen_filename, PersLen_dist,PersLen_cal_tors, append_dumps);
        Dumps.push_back(DLB);

        geninfile->add_entry(GENINFILE_DUMPS,"LBn"    ,PersLen_dump_every);
//        geninfile->add_entry(GENINFILE_DUMPS,"LBfn"   ,PersLen_filename);
        geninfile->add_entry(GENINFILE_DUMPS,"LBdist" ,PersLen_dist);
        geninfile->add_entry(GENINFILE_DUMPS,"LBtors" ,PersLen_cal_tors);
        geninfile->add_newline(GENINFILE_DUMPS);
    }

    /*
        Dump Plectoneme Statistics
    */
    int  PlecStats_dump_every      = 0;
    std::string PlecStats_filename = dump_dir+".plecstats";
    PlecStats_dump_every = InputChoice_get_single<int>         ("PlecStatsn"    ,&input,argv,PlecStats_dump_every);
    PlecStats_filename   = InputChoice_get_single<std::string> ("PlecStatsfn"   ,&input,argv,PlecStats_filename);

    if (PlecStats_dump_every>0) {
        Dump_PlecStats* DPlecStats = new Dump_PlecStats(_chain, PlecStats_dump_every, PlecStats_filename,append_dumps);
        Dumps.push_back(DPlecStats);

        geninfile->add_entry(GENINFILE_DUMPS,"PlecStatsn"  ,PlecStats_dump_every);
//        geninfile->add_entry(GENINFILE_DUMPS,"PlecStatsfn" ,PlecStats_filename);
        geninfile->add_newline(GENINFILE_DUMPS);
    }

    /*
        Dump ProximityPlec xp
    */
    int  xp_dump_every      = 0;
    std::string xp_filename = dump_dir+".xp";
    double  xp_threshold_dist = 0;
    int     xp_min_seg_dist   = 0;
    double  frac_closest      = 0.75;

    xp_dump_every = InputChoice_get_single<int>         ("XPn"    ,&input,argv,xp_dump_every);
    xp_dump_every = InputChoice_get_single<int>         ("xpn"    ,&input,argv,xp_dump_every);
    xp_filename   = InputChoice_get_single<std::string> ("XPfn"   ,&input,argv,xp_filename);
    xp_filename   = InputChoice_get_single<std::string> ("xpfn"   ,&input,argv,xp_filename);

    xp_threshold_dist = InputChoice_get_single<double> ("XPthresdist"  ,&input,argv,xp_threshold_dist);
    xp_min_seg_dist   = InputChoice_get_single<double> ("XPminsegdist" ,&input,argv,xp_min_seg_dist);
    frac_closest      = InputChoice_get_single<double> ("XPclfrac"     ,&input,argv,frac_closest);


    if (xp_dump_every>0 && xp_threshold_dist > 0 && xp_min_seg_dist > 0) {
        Dump_ProximityPlec* Dxp= new Dump_ProximityPlec(_chain, xp_dump_every, xp_filename,xp_threshold_dist,xp_min_seg_dist,frac_closest,append_dumps);
        Dumps.push_back(Dxp);

        geninfile->add_entry(GENINFILE_DUMPS,"XPn"  ,xp_dump_every);
//        geninfile->add_entry(GENINFILE_DUMPS,"PlecStatsfn" ,PlecStats_filename);
        geninfile->add_newline(GENINFILE_DUMPS);
    }

    /*
        Dump Restart Snapshots
    */
    int    Restart_dump_every    = 0;
    std::string Restart_filename = dump_dir+".restart";

    Restart_dump_every = InputChoice_get_single<int>         ("Restartn"    ,&input,argv,Restart_dump_every);
    Restart_dump_every = InputChoice_get_single<int>         ("restartn"    ,&input,argv,Restart_dump_every);
    Restart_filename   = InputChoice_get_single<std::string> ("Restartfn"   ,&input,argv,Restart_filename);
    Restart_filename   = InputChoice_get_single<std::string> ("restartfn"   ,&input,argv,Restart_filename);

    if (Restart_dump_every>0) {
        Dump_Restart* DRestart = new Dump_Restart(_chain, Restart_dump_every, Restart_filename,append_dumps);
        Dumps.push_back(DRestart);

        geninfile->add_entry(GENINFILE_DUMPS,"Restartn"  ,Restart_dump_every);
//        geninfile->add_entry(GENINFILE_DUMPS,"Restartfn" ,Restart_filename);
        geninfile->add_newline(GENINFILE_DUMPS);
    }

    /*
        Dump Solenoid Correlation Function
    */
    int    SolCor_dump_every        = 0;
    std::string SolCor_filename     = dump_dir+".solcor";
    int    SolCor_k_max             = 100;
    int    SolCor_n_max             = 1;

    SolCor_dump_every = InputChoice_get_single<int>         ("SCn"    ,&input,argv,SolCor_dump_every);
    SolCor_filename   = InputChoice_get_single<std::string> ("SCfn"   ,&input,argv,SolCor_filename);
    SolCor_k_max      = InputChoice_get_single<int>         ("SCkmax" ,&input,argv,SolCor_k_max);
    SolCor_n_max      = InputChoice_get_single<int>         ("SCnmax" ,&input,argv,SolCor_n_max);

    if (SolCor_dump_every>0) {
        arma::colvec z_dir = {0,0,1};
        Dump_SolCor* DSC = new Dump_SolCor(_chain, SolCor_dump_every, SolCor_filename, z_dir, SolCor_k_max, SolCor_n_max, append_dumps);
        Dumps.push_back(DSC);

        geninfile->add_entry(GENINFILE_DUMPS,"SCn"    ,SolCor_dump_every);
//        geninfile->add_entry(GENINFILE_DUMPS,"SCfn"   ,SolCor_filename);
        geninfile->add_entry(GENINFILE_DUMPS,"SCkmax" ,SolCor_k_max);
        geninfile->add_entry(GENINFILE_DUMPS,"SCnmax" ,SolCor_n_max);
        geninfile->add_newline(GENINFILE_DUMPS);
    }

    /*
        Dump State
    */
    int    State_dump_every    = 0;
    std::string State_filename = dump_dir+".state";
    bool   State_dump_triads   = false;
    bool   State_dump_Omegas   = false;

    State_dump_every = InputChoice_get_single<int>         ("Stn"    ,&input,argv,State_dump_every);
    State_filename   = InputChoice_get_single<std::string> ("Stfn"   ,&input,argv,State_filename);

    State_dump_triads          = parse_flag(State_dump_triads, "-Sttriads",argv);
    State_dump_triads          = parse_flag(State_dump_triads, "-Sttds"   ,argv);
    State_dump_Omegas          = parse_flag(State_dump_Omegas, "-StOmegas",argv);
    State_dump_Omegas          = parse_flag(State_dump_Omegas, "-StOm"    ,argv);

    if (!State_dump_triads) {
        State_dump_triads = InputChoice_get_single<bool> ("Sttriads" ,&input,argv,State_dump_triads);
        State_dump_triads = InputChoice_get_single<bool> ("Sttds"    ,&input,argv,State_dump_triads);
    }
    if (!State_dump_Omegas) {
        State_dump_Omegas = InputChoice_get_single<bool> ("StOmegas" ,&input,argv,State_dump_Omegas);
        State_dump_Omegas = InputChoice_get_single<bool> ("StOm"     ,&input,argv,State_dump_Omegas);
    }

    if (State_dump_every>0) {
        Dump_State* DState = new Dump_State(_chain, mode, State_dump_every, State_filename,State_dump_triads,State_dump_Omegas,EV_rad, append_dumps);
        Dumps.push_back(DState);

        geninfile->add_entry(GENINFILE_DUMPS,"Stn"      ,State_dump_every);
//        geninfile->add_entry(GENINFILE_DUMPS,"Stfn"     ,State_filename);
        geninfile->add_entry(GENINFILE_DUMPS,"Sttriads" ,State_dump_triads);
        geninfile->add_entry(GENINFILE_DUMPS,"StOmegas" ,State_dump_Omegas);
        geninfile->add_newline(GENINFILE_DUMPS);
    }

    /*
        Dump Stiff
    */
    int    Stiff_dump_every = 0;
    std::string Stiff_filename   = dump_dir+".stiff";
    bool   Stiff_all_bps    = false;

    Stiff_dump_every = InputChoice_get_single<int>         ("Stiffn"    ,&input,argv,Stiff_dump_every);
    Stiff_filename   = InputChoice_get_single<std::string> ("Stifffn"   ,&input,argv,Stiff_filename);

    Stiff_all_bps    = parse_flag(false, "-stiffall",argv);
    Stiff_all_bps    = parse_flag(false, "-Stiffall",argv);
    if (!Stiff_all_bps) {
        Stiff_all_bps = InputChoice_get_single<bool>         ("Stiffall" ,&input,argv,Stiff_all_bps);
    }

    if (Stiff_dump_every>0) {
        Dump_Stiff* DStiff = new Dump_Stiff(_chain, Stiff_dump_every, Stiff_filename, Stiff_all_bps);
        Dumps.push_back(DStiff);

        geninfile->add_entry(GENINFILE_DUMPS,"Stiffn"   ,Stiff_dump_every);
//        geninfile->add_entry(GENINFILE_DUMPS,"Stifffn"  ,Stiff_filename);
        geninfile->add_entry(GENINFILE_DUMPS,"Stiffall" ,Stiff_all_bps);
        geninfile->add_newline(GENINFILE_DUMPS);
    }

    /*
        Dump Supercoiling Phase
    */
    int    SupPhase_dump_every = 0;
    std::string SupPhase_filename   = dump_dir+".scp";
    SupPhase_dump_every = InputChoice_get_single<int>         ("SCPn"    ,&input,argv,SupPhase_dump_every);
    SupPhase_filename   = InputChoice_get_single<std::string> ("SCPfn"   ,&input,argv,SupPhase_filename);
    if (SupPhase_dump_every>0) {
        Dump_SupercoilingPhase* DSupPhase = new Dump_SupercoilingPhase(_chain, SupPhase_dump_every, SupPhase_filename ,append_dumps);
        Dumps.push_back(DSupPhase);

        geninfile->add_entry(GENINFILE_DUMPS,"SCPn"  ,SupPhase_dump_every);
//        geninfile->add_entry(GENINFILE_DUMPS,"SCPfn" ,SupPhase_filename);
        geninfile->add_newline(GENINFILE_DUMPS);
    }

    /*
        Dump Tangent Correlation
    */
    int    TanCor_dump_every = 0;
    std::string TanCor_filename   = dump_dir+".tancor";
    int    TanCor_dist       = 200;

    TanCor_dump_every = InputChoice_get_single<int>         ("TCn"    ,&input,argv,TanCor_dump_every);
    TanCor_filename   = InputChoice_get_single<std::string> ("TCfn"   ,&input,argv,TanCor_filename);
    TanCor_dist       = InputChoice_get_single<int>         ("TCdist" ,&input,argv,TanCor_dist);

    if (TanCor_dump_every>0) {
        Dump_TanCor* DTC = new Dump_TanCor(_chain, TanCor_dump_every, TanCor_filename, TanCor_dist, append_dumps);
        Dumps.push_back(DTC);

        geninfile->add_entry(GENINFILE_DUMPS,"TCn"    ,TanCor_dump_every);
//        geninfile->add_entry(GENINFILE_DUMPS,"TCfn"   ,TanCor_filename);
        geninfile->add_entry(GENINFILE_DUMPS,"TCdist" ,TanCor_dist);
        geninfile->add_newline(GENINFILE_DUMPS);
    }

    /*
        Dump Thetas
    */
    int    Thetas_dump_every = 0;
    std::string Thetas_filename   = dump_dir+".thetas";

    Thetas_dump_every = InputChoice_get_single<int>         ("Thetasn"    ,&input,argv,Thetas_dump_every);
    Thetas_filename   = InputChoice_get_single<std::string> ("Thetasfn"   ,&input,argv,Thetas_filename);

    if (Thetas_dump_every>0) {
        Dump_Thetas* DThetas = new Dump_Thetas(_chain, Thetas_dump_every, Thetas_filename ,append_dumps);
        Dumps.push_back(DThetas);

        geninfile->add_entry(GENINFILE_DUMPS,"Thetasn"  ,Thetas_dump_every);
//        geninfile->add_entry(GENINFILE_DUMPS,"Thetasfn" ,Thetas_filename);
        geninfile->add_newline(GENINFILE_DUMPS);
    }

    /*
        Dump AvgThetas
    */
    int    AvgThetas_dump_every = 0;
    std::string AvgThetas_filename   = dump_dir+".avgthetas";

    AvgThetas_dump_every = InputChoice_get_single<int>         ("AvgThetasn"    ,&input,argv,AvgThetas_dump_every);
    AvgThetas_filename   = InputChoice_get_single<std::string> ("AvgThetasfn"   ,&input,argv,AvgThetas_filename);

    if (AvgThetas_dump_every>0) {
        Dump_avgThetas* DAvgThetas = new Dump_avgThetas(_chain, AvgThetas_dump_every, AvgThetas_filename ,append_dumps);
        Dumps.push_back(DAvgThetas);

        geninfile->add_entry(GENINFILE_DUMPS,"AvgThetasn"  ,AvgThetas_dump_every);
//        geninfile->add_entry(GENINFILE_DUMPS,"AvgThetasfn" ,AvgThetas_filename);
        geninfile->add_newline(GENINFILE_DUMPS);
    }
    /*
        Dump Torque
    */
    int    Torque_dump_every = 0;
    std::string Torque_filename   = dump_dir+".torque";

    Torque_dump_every = InputChoice_get_single<int>         ("TRQn"    ,&input,argv,Torque_dump_every);
    Torque_filename   = InputChoice_get_single<std::string> ("TRQfn"   ,&input,argv,Torque_filename);

    if (Torque_dump_every>0) {
        Dump_Torque* DTRQ = new Dump_Torque(_chain, Torque_dump_every, Torque_filename ,append_dumps);
        Dumps.push_back(DTRQ);

        geninfile->add_entry(GENINFILE_DUMPS,"TRQn"  ,Torque_dump_every);
//        geninfile->add_entry(GENINFILE_DUMPS,"TRQfn" ,Torque_filename);
        geninfile->add_newline(GENINFILE_DUMPS);
    }


    /*
        Dump TorqueTrap
    */
    int    TorqueTrap_dump_every = 0;
    std::string TorqueTrap_filename   = dump_dir+".torquetrap";

    TorqueTrap_dump_every = InputChoice_get_single<int>         ("TRQTRPn"    ,&input,argv,TorqueTrap_dump_every);
    TorqueTrap_filename   = InputChoice_get_single<std::string> ("TRQTRPfn"   ,&input,argv,TorqueTrap_filename);

    if (TorqueTrap_dump_every>0) {
        Dump_TorqueTrap* DTRQTRP = new Dump_TorqueTrap(_chain, TorqueTrap_dump_every, TorqueTrap_filename ,append_dumps);
        Dumps.push_back(DTRQTRP);

        geninfile->add_entry(GENINFILE_DUMPS,"TRQTRPn"  ,TorqueTrap_dump_every);
//        geninfile->add_entry(GENINFILE_DUMPS,"TRQTRPfn" ,TorqueTrap_filename);
        geninfile->add_newline(GENINFILE_DUMPS);
    }

    /*
        Dump Torsional Stiffness
    */
    int    Ceff_dump_every    = 0;
    double Ceff_frac          = 1;
    std::string Ceff_filename      = dump_dir+".ceff";
    std::string Ceff_writhe_method = "fuller";

    Ceff_dump_every    = InputChoice_get_single<int>         ("Ceffn"    ,&input,argv,Ceff_dump_every);
    Ceff_dump_every    = InputChoice_get_single<int>         ("Ceffn"    ,&input,argv,Ceff_dump_every);
    Ceff_frac          = InputChoice_get_single<double>      ("Cefffrac" ,&input,argv,Ceff_frac);
    Ceff_frac          = InputChoice_get_single<double>      ("cefffrac" ,&input,argv,Ceff_frac);
    Ceff_filename      = InputChoice_get_single<std::string> ("Cefffn"   ,&input,argv,Ceff_filename);
    Ceff_filename      = InputChoice_get_single<std::string> ("cefffn"   ,&input,argv,Ceff_filename);
    Ceff_writhe_method = InputChoice_get_single<std::string> ("CeffWr"   ,&input,argv,Ceff_writhe_method);
    Ceff_writhe_method = InputChoice_get_single<std::string> ("Ceffwr"   ,&input,argv,Ceff_writhe_method);
    Ceff_writhe_method = InputChoice_get_single<std::string> ("ceffwr"   ,&input,argv,Ceff_writhe_method);


    if (Ceff_dump_every>0) {
        if (Ceff_frac > 1) Ceff_frac = 1;
        if (Ceff_frac < 0) Ceff_frac = 0;
        int start_id  = (1-Ceff_frac)*0.5*num_bp;
        int end_id    = (0.5 + Ceff_frac*0.5) *num_bp;

        Dump_EffTorsStiff* DCeff = new Dump_EffTorsStiff(_chain, Ceff_dump_every, Ceff_filename,_chain->get_force_dir(), Ceff_writhe_method, start_id, end_id ,append_dumps);
        Dumps.push_back(DCeff);

        geninfile->add_entry(GENINFILE_DUMPS,"Ceffn"    ,Ceff_dump_every);
        geninfile->add_entry(GENINFILE_DUMPS,"Cefffrac" ,Ceff_frac);
//        geninfile->add_entry(GENINFILE_DUMPS,"Cefffn"   ,Ceff_filename);
        geninfile->add_entry(GENINFILE_DUMPS,"CeffWr"   ,Ceff_writhe_method);
        geninfile->add_newline(GENINFILE_DUMPS);
    }

    /*
        Dump Twist
    */
    int         TW_dump_every = 0;
    std::string TW_filename   = dump_dir+".tw";
    TW_dump_every  = InputChoice_get_single<int>         ("Twn" ,&input,argv,TW_dump_every);
    TW_dump_every  = InputChoice_get_single<int>         ("twn" ,&input,argv,TW_dump_every);
    TW_dump_every  = InputChoice_get_single<int>         ("TWn" ,&input,argv,TW_dump_every);
    TW_filename    = InputChoice_get_single<std::string> ("Twfn",&input,argv,TW_filename);
    TW_filename    = InputChoice_get_single<std::string> ("twfn",&input,argv,TW_filename);
    TW_filename    = InputChoice_get_single<std::string> ("TWfn",&input,argv,TW_filename);

    if (TW_dump_every>0) {
        Dump_Twist* DTW = new Dump_Twist(_chain, TW_dump_every, TW_filename,append_dumps);
        Dumps.push_back(DTW);

        geninfile->add_entry(GENINFILE_DUMPS,"Twn" ,TW_dump_every);
        geninfile->add_newline(GENINFILE_DUMPS);
    }

    /*
        Dump Writhe Map
    */
    int    WM_dump_every = 0;
    std::string WM_filename   = dump_dir+".wm";
    int    WM_seg_size   = 1;

    WM_dump_every    = InputChoice_get_single<int>         ("WMn"    ,&input,argv,WM_dump_every);
    WM_filename      = InputChoice_get_single<std::string> ("WMfn"   ,&input,argv,WM_filename);
    WM_seg_size      = InputChoice_get_single<int>         ("WMseg"  ,&input,argv,WM_seg_size);

    if (WM_dump_every>0) {
        if (WM_seg_size < 1) {
            WM_seg_size = 1;
        }
        if (WM_seg_size > num_bp/10) {
            WM_seg_size = num_bp/10;
        }
        Dump_WritheMap* DWM = new Dump_WritheMap(_chain, WM_dump_every, WM_filename, WM_seg_size, append_dumps);
        Dumps.push_back(DWM);

        geninfile->add_entry(GENINFILE_DUMPS,"WMn"   ,WM_dump_every);
//        geninfile->add_entry(GENINFILE_DUMPS,"WMfn"  ,WM_filename);
        geninfile->add_entry(GENINFILE_DUMPS,"WMseg" ,WM_seg_size);
        geninfile->add_newline(GENINFILE_DUMPS);
    }

    /*
        Dump XYZ
    */
    int    XYZ_dump_every      = 0;
    std::string XYZ_dumpxyzfilename = dump_dir+".xyz";
    std::string XYZ_dumpxyzcenter   = "first";
    std::string XYZ_rep             = "simple"; //"EV"

    XYZ_dump_every      = InputChoice_get_single<int>         ("XYZn"          ,&input,argv,XYZ_dump_every);
    XYZ_dumpxyzfilename = InputChoice_get_single<std::string> ("XYZfn"         ,&input,argv,XYZ_dumpxyzfilename);
    XYZ_dumpxyzcenter   = InputChoice_get_single<std::string> ("XYZ_translate" ,&input,argv,XYZ_dumpxyzcenter);
    XYZ_rep             = InputChoice_get_single<std::string> ("XYZ_repr"      ,&input,argv,XYZ_rep);

    if (XYZ_dump_every>0) {
        Dump_xyz* Dxyz = new Dump_xyz(_chain, XYZ_dump_every, XYZ_dumpxyzfilename ,append_dumps, XYZ_dumpxyzcenter,XYZ_rep);
        Dumps.push_back(Dxyz);

        geninfile->add_entry(GENINFILE_DUMPS,"XYZn"          ,XYZ_dump_every);
//        geninfile->add_entry(GENINFILE_DUMPS,"XYZfn"         ,XYZ_dumpxyzfilename);
        geninfile->add_entry(GENINFILE_DUMPS,"XYZ_translate" ,XYZ_dumpxyzcenter);
        geninfile->add_entry(GENINFILE_DUMPS,"XYZ_repr"      ,XYZ_rep);
        geninfile->add_newline(GENINFILE_DUMPS);
    }

    /*
        Dump XYZ extreme extension
    */
    int         XYZ_extreme_ext_dump_every      = 0;
    std::string XYZ_extreme_ext_dumpxyzfilename = dump_dir+"_extreme_ext";
    std::string XYZ_extreme_ext_dumpxyzcenter   = "first";
    std::string XYZ_extreme_ext_rep             = "simple"; //"EV"

    int         XYZ_extreme_ext_min_datapoints = 1000000;
    double      XYZ_extreme_ext_num_sigmas     = 3;

    XYZ_extreme_ext_dump_every      = InputChoice_get_single<int>         ("XYZxtrzn"          ,&input,argv,XYZ_extreme_ext_dump_every);
    XYZ_extreme_ext_dumpxyzfilename = InputChoice_get_single<std::string> ("XYZxtrzfn"         ,&input,argv,XYZ_extreme_ext_dumpxyzfilename);
    XYZ_extreme_ext_dumpxyzcenter   = InputChoice_get_single<std::string> ("XYZxtrz_translate" ,&input,argv,XYZ_extreme_ext_dumpxyzcenter);
    XYZ_extreme_ext_rep             = InputChoice_get_single<std::string> ("XYZxtrz_repr"      ,&input,argv,XYZ_extreme_ext_rep);

    XYZ_extreme_ext_min_datapoints  = InputChoice_get_single<int>         ("XYZxtrz_min_datapoints",&input,argv,XYZ_extreme_ext_min_datapoints);
    XYZ_extreme_ext_num_sigmas      = InputChoice_get_single<double>      ("XYZxtrz_num_sigma"     ,&input,argv,XYZ_extreme_ext_num_sigmas);

    if (XYZ_extreme_ext_dump_every>0) {
        Dump_xyz_extreme_ext* Dxyz_xtrz = new Dump_xyz_extreme_ext(_chain, XYZ_extreme_ext_dump_every, XYZ_extreme_ext_dumpxyzfilename, XYZ_extreme_ext_min_datapoints, XYZ_extreme_ext_num_sigmas,append_dumps, XYZ_extreme_ext_dumpxyzcenter,XYZ_extreme_ext_rep);
        Dumps.push_back(Dxyz_xtrz);

        geninfile->add_entry(GENINFILE_DUMPS,"XYZxtrzn"          ,XYZ_extreme_ext_dump_every);
//        geninfile->add_entry(GENINFILE_DUMPS,"XYZfn"         ,XYZ_dumpxyzfilename);
        geninfile->add_entry(GENINFILE_DUMPS,"XYZxtrz_translate" ,XYZ_extreme_ext_dumpxyzcenter);
        geninfile->add_entry(GENINFILE_DUMPS,"XYZxtrz_repr"      ,XYZ_extreme_ext_rep);
        geninfile->add_entry(GENINFILE_DUMPS,"XYZxtrz_min_datapoints" ,XYZ_extreme_ext_min_datapoints);
        geninfile->add_entry(GENINFILE_DUMPS,"XYZxtrz_num_sigma"      ,XYZ_extreme_ext_num_sigmas);
        geninfile->add_newline(GENINFILE_DUMPS);
    }

    /*
        Dump Triads
    */
    int    Triads_dump_every            = 0;
    std::string Triads_dumpxyzfilename  = dump_dir+".triads";
    std::string Triads_dumpxyzcenter    = "first";

    Triads_dump_every      = InputChoice_get_single<int>         ("Trdsn"          ,&input,argv,Triads_dump_every);
    Triads_dumpxyzfilename = InputChoice_get_single<std::string> ("Trdsfn"         ,&input,argv,Triads_dumpxyzfilename);
    Triads_dumpxyzcenter   = InputChoice_get_single<std::string> ("Trds_translate" ,&input,argv,Triads_dumpxyzcenter);

    if (Triads_dump_every>0) {
        Dump_xyz* Dxyz = new Dump_xyz(_chain, Triads_dump_every, Triads_dumpxyzfilename ,append_dumps, Triads_dumpxyzcenter,"pdb");
        Dumps.push_back(Dxyz);

        geninfile->add_entry(GENINFILE_DUMPS,"Trdsn"          ,Triads_dump_every);
//        geninfile->add_entry(GENINFILE_DUMPS,"Trdsfn"         ,Triads_dumpxyzfilename);
        geninfile->add_entry(GENINFILE_DUMPS,"Trds_translate" ,Triads_dumpxyzcenter);
        geninfile->add_newline(GENINFILE_DUMPS);
    }

    return Dumps;

}
