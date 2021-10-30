#ifndef __MCS_UNWINDER_H__
#define __MCS_UNWINDER_H__
#include "MCStep.h"
#include "MCS_CSrot.h"

#define MCSUW_CHECK_LINK_EVERY    10000
#define MCSUW_CHECK_LINK_SEG_SIZE 5

#define MCSUW_KILL_ON_DLK_VIOLATION 1


class MCS_Unwinder;

class MCS_Unwinder: public MCStep
{
protected:
    /*
        Winding info
    */
    int     steps_unwinding;
    int     steps_winding;
    int     steps_winding_per_split;
    int     step_equi;
    double  sigma_unwound;
    double  dLK_wound;
    double  dLK_unwound;
    int     n_split_winding;
    std::vector<double> winding_dLK;

    /*
        MC Steps
    */
    std::vector<MCStep*> MCSteps;

    /*
        Backup Configurations
    */
    arma::mat   backup_bp_pos;
    arma::cube  backup_triads;

    double link_check_density;


public:
    MCS_Unwinder(Chain * ch,const std::vector<long long int> & seedseq, double sigma_unwound, int steps_unwinding, int steps_winding, int step_equi, int n_split_winding);
    ~MCS_Unwinder();

    void set_excluded_volume(ExVol* EV);

    void update_settings();

// Monte Carlo Step (call main functionality)
public:
    bool MC();
    bool MC_move();

    void run(long long int steps);


};

#endif
