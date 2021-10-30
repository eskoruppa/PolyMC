#ifndef __MCS_REP_THERM_H__
#define __MCS_REP_THERM_H__
#include "MCStep.h"
#include "MCS_CSrot.h"

#define LEVER_SIZE_FAC  10 //
#define MIN_INTSEG_SIZE 2
#define CAP_RETRIES 1000



#define SELECT_EXTRA_CRITERION_FACTOR       1

#define REM_NEWTON_RAPHSON_MAX_STEPS        25
#define REM_NEWTON_RAPHSON_EPS              1e-15
#define REM_DISPLACEMENT_EPS                1e-10

#define ADD_NEWTON_RAPHSON_THETA_EPS        1e-15
#define ADD_NEWTON_RAPHSON_DISPLACEMENT_EPS 1e-10
#define ADD_NEWTON_RAPHSON_MAX_STEPS        25

#define THERMALIZE_DIST                     6
#define THERMALIZE_MOVES_PER_RESIDUE        4

#define THERMALIZE_S_MIN                    2



/*
ToDo: Build in a check for the Newton-Raphson elongator to identify any of the tangents has
an orthogonal component that is oppositvely oriented to the vector w. In that case This
method should not be used and rather a simple pivot rotatation should be executed.

*/




class MCS_RepTherm;

class MCS_RepTherm: public MCStep
{
protected:
    bool   closed;

    MCS_CSrot * MC_pre_thermal_rot;
    MCS_CSrot * MC_thermal_rot;

    int rep_size_min,rep_size_max,rep_dist_min,rep_dist_max;
    std::uniform_int_distribution<> gen_idFirst{0, 1};
    std::uniform_int_distribution<> gen_repsize{0, 1};
    std::uniform_int_distribution<> gen_repdist{0, 1};

    bool twist_translation_active;
    bool thermalize;

    arma::mat   backup_bp_pos1;
    arma::cube  backup_triads1;
    arma::mat   backup_bp_pos2;
    arma::cube  backup_triads2;

    int idFirst,idLast;
    int idThermFirst,idThermLast;
    int idEnergyFirst,idEnergyLast;

    int rem_hA,rem_hB,rem_hAp,rem_hBp;
    int rem_hv1,rem_hv2,rem_hvp;
    int add_hA,add_hB,add_hAp,add_hBp;
    int add_hv,add_hv1p,add_hv2p;
    int rep_size,rep_dist;
    int lever_size,intseg_size;
    int total_size;

    bool rep_forward;

//    arma::colvec rem_v1,rem_v2, rem_w;
//    double len_rem_v1,len_rem_v2,len_rem_w;

//    arma::colvec add_v1,add_v2, add_w;
//    double len_add_v1,len_add_v2,len_add_w;


    /*
        This keeps track the the average theta0 which will be used as a starting
        value for the Newton Raphson minimization.
    */
    double add_NR_mean_theta0=M_PI*(0.22);
    int    add_NR_count_theta0=1;

    double avg_deltaE    = 0;
    int avg_deltaE_count = 0;



public:
    MCS_RepTherm( Chain * ch, int rep_size_min, int rep_size_max, int rep_dist_min, int rep_dist_max,bool twist_translation_active=true, bool thermalize=true);
    ~MCS_RepTherm();

    void set_excluded_volume(ExVol* EVol);

    void update_settings();

// Monte Carlo Step (call main functionality)
public:
    bool MC();
    bool MC_move();


protected:
    void conf_make_backup(int idThermFirst, int idThermLast);
    void conf_revert_to_backup(int idThermFirst, int idThermLast);

    void translate_twist();
    bool REM_move();
    void MID_move();
    bool ADD_move();
    bool ADD_move_FullNewtonRaphson();

    bool Pre_Thermalize();
    bool Thermalize();
    bool Thermalize_Hinges();
    bool Check_EV();



};

#endif
