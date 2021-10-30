#ifndef __PAIR_LJ_INCLUDED__
#define __PAIR_LJ_INCLUDED__

#include "../Pair.h"


class Pair_LJ;

class Pair_LJ : public Pair {

protected:
    double epsilon;
    double rmin;

public:

    Pair_LJ(Chain * ch,
            TypeGroup * type_group_A,
            TypeGroup * type_group_B,
            double cutoff,
            std::vector<double> params);
    Pair_LJ();
    virtual ~Pair_LJ() override;

    /*---------------------------------------------------------------------------------------------------------------------*/

    virtual double Eval_change_type(const arma::mat * pos) override;
    virtual double Eval_change_type(const arma::mat * pos,
                            int monomer_id,
                            int type_from,
                            int type_to) override;

    /*---------------------------------------------------------------------------------------------------------------------*/

    virtual double Eval_Energy(const arma::mat * pos) override;

    virtual double Eval_Delta_Energy(   const arma::mat * pos_new,
                                        const arma::mat * pos_old) override;
    virtual double Eval_Delta_Energy(   const arma::mat * pos_new,
                                        const arma::mat * pos_old,
                                        const std::vector<int> * A_unmoved,
                                        const std::vector<int> * A_moved,
                                        const std::vector<int> * A_individual,
                                        const std::vector<int> * B_unmoved,
                                        const std::vector<int> * B_moved,
                                        const std::vector<int> * B_individual) override;

    /*---------------------------------------------------------------------------------------------------------------------*/


protected:

    double Potential(double distance);

//    int find_id_in_group_A(int monomer_id);
//    int find_id_in_group_B(int monomer_id);

};





#endif
