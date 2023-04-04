#include "Pair.h"

#include "../SO3Methods.h"
#include "../ExtraFuncs.h"


Pair::Pair( Chain * ch,
            TypeGroup * type_group_A,
            TypeGroup * type_group_B,
            double cutoff,
            std::vector<double> params) :
chain(ch),
num_bp(ch->get_num_bp()),
disc_len(ch->get_disc_len()),
type_group_A(type_group_A),
type_group_B(type_group_B),
type_A(&type_group_A->type),
type_B(&type_group_B->type),
type_id_A(type_A->id),
type_id_B(type_B->id),
group_A(&type_group_A->elements),
group_B(&type_group_B->elements),
params(params),
A_unmoved(&type_group_A->unmoved),
A_moved(&type_group_A->moved),
A_individual(&type_group_A->individual),
B_unmoved(&type_group_B->unmoved),
B_moved(&type_group_B->moved),
B_individual(&type_group_B->individual),
num_moved_A(&type_group_A->num_moved),
num_moved_B(&type_group_B->num_moved)
{
    if (type_id_A == type_id_B) single_group = true;
    else                        single_group = false;

    upper_cutoff = cutoff;
    lower_cutoff = type_A->hard_repulsion_cuttoff + type_B->hard_repulsion_cuttoff;

    type_group_A->involved_in_pairinteraction = true;
}




Pair::Pair() {

}

Pair::~Pair() {

}


/*---------------------------------------------------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------------------------------------------------*/


double Pair::Eval_change_type(const arma::mat * pos) {

    double Delta_E = 0;
    double dist;
    arma::colvec new_pos;

    if (single_group) {

        /* ------------------------------ */
        /*          OLD ENERGY            */
        /* ------------------------------ */
        for (unsigned i=0;i<type_group_A->elements_removed.size();i++) {
            // removed A with remaining A
            for (unsigned j=0;j<type_group_A->elements_remaining.size();j++) {
                dist = arma::norm( pos->col(type_group_A->elements_removed[i]) - pos->col(type_group_A->elements_remaining[j]) );
                #ifdef PAIR_DEBUG
                if (dist < lower_cutoff) {
                    std::cout << "Error: Pair detected violation of hard core repulsion!" << std::endl;
                    std::exit(0);
                }
                #endif
                if (dist < upper_cutoff) {
                    Delta_E -= Potential(dist);
                }
            }
            // removed A with removed A
            for (unsigned j=i+1;j<type_group_A->elements_removed.size();j++) {
                dist = arma::norm( pos->col(type_group_A->elements_removed[i]) - pos->col(type_group_A->elements_removed[j]) );
                #ifdef PAIR_DEBUG
                if (dist < lower_cutoff) {
                    std::cout << "Error: Pair detected violation of hard core repulsion!" << std::endl;
                    std::exit(0);
                }
                #endif
                if (dist < upper_cutoff) {
                    Delta_E -= Potential(dist);
                }
            }
        }

        /* ------------------------------ */
        /*          NEW ENERGY            */
        /* ------------------------------ */
        for (unsigned i=0;i<type_group_A->elements_added.size();i++) {
            // added A with remaining A
            for (unsigned j=0;j<type_group_A->elements_remaining.size();j++) {
                dist = arma::norm( pos->col(type_group_A->elements_added[i]) - pos->col(type_group_A->elements_remaining[j]) );
                if (dist < lower_cutoff) {
                    return PAIR_INFINITY;
                }
                else if (dist < upper_cutoff) {
                    Delta_E += Potential(dist);
                }
            }
             // added A with added A
            for (unsigned j=i+1;j<type_group_A->elements_added.size();j++) {
                dist = arma::norm( pos->col(type_group_A->elements_added[i]) - pos->col(type_group_A->elements_added[j]) );
                if (dist < lower_cutoff) {
                    return PAIR_INFINITY;
                }
                else if (dist < upper_cutoff) {
                    Delta_E += Potential(dist);
                }
            }
        }
    }
    else {

        /* ------------------------------ */
        /*          OLD ENERGY            */
        /* ------------------------------ */

        // removed A with old B
        for (unsigned i=0;i<type_group_A->elements_removed.size();i++) {
            for (unsigned j=0;j<type_group_B->elements.size();j++) {
                dist = arma::norm( pos->col(type_group_A->elements_removed[i]) - pos->col(type_group_B->elements[j]) );
                #ifdef PAIR_DEBUG
                if (dist < lower_cutoff) {
                    std::cout << "Error: Pair detected violation of hard core repulsion!" << std::endl;
                    std::exit(0);
                }
                #endif
                if (dist < upper_cutoff) {
                    Delta_E -= Potential(dist);
                }
            }
        }

        // removed B with old A
        for (unsigned i=0;i<type_group_B->elements_removed.size();i++) {
            for (unsigned j=0;j<type_group_A->elements.size();j++) {
                dist = arma::norm( pos->col(type_group_B->elements_removed[i]) - pos->col(type_group_A->elements[j]) );
                #ifdef PAIR_DEBUG
                if (dist < lower_cutoff) {
                    std::cout << "Error: Pair detected violation of hard core repulsion!" << std::endl;
                    std::exit(0);
                }
                #endif
                if (dist < upper_cutoff) {
                    Delta_E -= Potential(dist);
                }
            }
        }

        /* ------------------------------ */
        /*          NEW ENERGY            */
        /* ------------------------------ */

        // added A with trial B
        for (unsigned i=0;i<type_group_A->elements_added.size();i++) {
            for (unsigned j=0;j<type_group_B->elements_trial.size();j++) {
                dist = arma::norm( pos->col(type_group_A->elements_added[i]) - pos->col(type_group_B->elements_trial[j]) );
                if (dist < lower_cutoff) {
                    return PAIR_INFINITY;
                }
                else if (dist < upper_cutoff) {
                    Delta_E += Potential(dist);
                }
            }
        }
        // added B with trial A
        for (unsigned i=0;i<type_group_B->elements_added.size();i++) {
            for (unsigned j=0;j<type_group_A->elements_trial.size();j++) {
                dist = arma::norm( pos->col(type_group_B->elements_added[i]) - pos->col(type_group_A->elements_trial[j]) );
                if (dist < lower_cutoff) {
                    return PAIR_INFINITY;
                }
                else if (dist < upper_cutoff) {
                    Delta_E += Potential(dist);
                }
            }
        }
    }
    return Delta_E;
}


double Pair::Eval_change_type(  const arma::mat * pos,
                                int monomer_id,
                                int type_from,
                                int type_to)
{

    #ifdef PAIR_DEBUG
    if ((type_from == type_id_A && find_id_in_group_A(monomer_id) == -1) || (type_from == type_id_B && find_id_in_group_B(monomer_id) == -1)) {
        std::cout << "Error: Pair detected invalid type change!" << std::endl;
        std::exit(0);
    }
    #endif

    double Delta_E = 0;
    double dist;
    arma::colvec new_pos;

    /* ------------------------------ */
    /*          OLD ENERGY            */
    /* ------------------------------ */

    if (type_from == type_id_A) {
        new_pos = pos->col(monomer_id);
        for (int i=0;i<group_B->size();i++) {
            dist = arma::norm( new_pos - pos->col(group_B->at(i)) );

            #ifdef PAIR_DEBUG
            if (dist < lower_cutoff) {
                std::cout << "Error: Pair detected violation of hard core repulsion!" << std::endl;
                std::exit(0);
            }
            #endif

            if (dist < upper_cutoff) {
                Delta_E -= Potential(dist);
            }
        }
    }
    else if (type_from == type_id_B) {
        new_pos = pos->col(monomer_id);
        for (int i=0;i<group_A->size();i++) {
            dist = arma::norm( new_pos - pos->col(group_A->at(i)) );

            #ifdef PAIR_DEBUG
            if (dist < lower_cutoff) {
                std::cout << "Error: Pair detected violation of hard core repulsion!" << std::endl;
                std::exit(0);
            }
            #endif

            if (dist < upper_cutoff) {
                Delta_E -= Potential(dist);
            }
        }
    }

    /* ------------------------------ */
    /*          NEW ENERGY            */
    /* ------------------------------ */

    if (type_to == type_id_A) {
        new_pos = pos->col(monomer_id);
        for (int i=0;i<group_B->size();i++) {
            dist = arma::norm( new_pos - pos->col(group_B->at(i)) );
            if (dist < lower_cutoff) {
                return PAIR_INFINITY;
            }
            else if (dist < upper_cutoff) {
                Delta_E += Potential(dist);
            }
        }
    }

    if (type_to == type_id_B) {
        new_pos = pos->col(monomer_id);
        for (int i=0;i<group_A->size();i++) {
            dist = arma::norm( new_pos - pos->col(group_A->at(i)) );
            if (dist < lower_cutoff) {
                return PAIR_INFINITY;
            }
            else if (dist < upper_cutoff) {
                Delta_E += Potential(dist);
            }
        }
    }
    return Delta_E;
}

/*---------------------------------------------------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------------------------------------------------*/

double Pair::Eval_Energy(const arma::mat * pos)
{
    double E,dist;
    arma::colvec pos_i;

    E = 0;
    if (single_group) {
        for (int i=0;i<group_A->size();i++) {
            pos_i = pos->col(i);
            for (int j=i+1;j<group_B->size();j++) {
                dist = arma::norm( pos_i - pos->col(group_B->at(j)) );
                if (dist < lower_cutoff) {
                    std::cout << "dist = " << dist << std::endl;
                    return PAIR_INFINITY;
                }
                else if (dist < upper_cutoff) {
                    E += Potential(dist);
                }
            }
        }
    }
    else {
        for (int i=0;i<group_A->size();i++) {
            pos_i = pos->col(i);
            for (int j=0;j<group_B->size();j++) {
                dist = arma::norm( pos_i - pos->col(group_B->at(j)) );
                if (dist < lower_cutoff) {
                    std::cout << "dist = " << dist << std::endl;
                    std::cout << "Ai = " << group_A->at(i) << std::endl;
                    std::cout << "Bj = " << group_B->at(j) << std::endl;
                    return PAIR_INFINITY;
                }
                else if (dist < upper_cutoff) {
                    E += Potential(dist);
                }
            }
        }
    }
    return E;
}


double Pair::Eval_Delta_Energy( const arma::mat * pos_new,
                                const arma::mat * pos_old)
{
    if (&num_moved_A == 0 && &num_moved_B == 0) {
        return 0;
    }
    return Eval_Delta_Energy(   pos_new,
                                pos_old,
                                this->A_unmoved,
                                this->A_moved,
                                this->A_individual,
                                this->B_unmoved,
                                this->B_moved,
                                this->B_individual);
}


double Pair::Eval_Delta_Energy( const arma::mat * pos_new,
                                const arma::mat * pos_old,
                                const std::vector<int> * A_unmoved,
                                const std::vector<int> * A_moved,
                                const std::vector<int> * A_individual,
                                const std::vector<int> * B_unmoved,
                                const std::vector<int> * B_moved,
                                const std::vector<int> * B_individual)
{

    arma::colvec p_new, p_old;
    int num_unmov, num_moved, num_indiv;
    int A1,A2,B1,B2;
    int A1p,A2p,B1p,B2p;

    int ref1,ref2;

    int interv1_first;
    int interv1_last;

    double Delta_E = 0;
    double dist_old,dist_new;


    num_unmov = A_unmoved->size()/2;
    num_moved = A_moved->size()/2;
    num_indiv = A_individual->size()/2;

    /*
        Inter-interval checks
    */
    for (int mov=0;mov<num_moved;mov++){
        ref1 = mov*2;
        ref2 = ref1+1;
        A1 = A_moved->at(ref1);
        A2 = A_moved->at(ref2);
        B1 = B_moved->at(ref1);
        B2 = B_moved->at(ref2);

        /*
            Moved - Unmoved checks
        */
        for (int unm=0;unm<num_unmov;unm++){
            ref1 = unm*2;
            ref2 = ref1+1;
            A1p = A_unmoved->at(ref1);
            A2p = A_unmoved->at(ref2);
            B1p = B_unmoved->at(ref1);
            B2p = B_unmoved->at(ref2);

            // check A with Bp
            for (int A=A1;A<=A2;A++) {
                for (int Bp=B1p;Bp<=B2p;Bp++) {
                    dist_new = arma::norm( pos_new->col(group_A->at(A)) - pos_new->col(group_B->at(Bp)) );
                    if (dist_new < lower_cutoff) {
                        return PAIR_INFINITY;
                    }
                    else if (dist_new < upper_cutoff) {
                        Delta_E += Potential(dist_new);
                    }

                    dist_old = arma::norm( pos_old->col(group_A->at(A)) - pos_old->col(group_B->at(Bp)) );
                    if (dist_old < upper_cutoff) {
                        Delta_E -= Potential(dist_old);
                    }
                }
            }

            // check B with Ap
            for (int B=B1;B<=B2;B++) {
                for (int Ap=A1p;Ap<=A2p;Ap++) {
                    dist_new = arma::norm( pos_new->col(group_A->at(Ap)) - pos_new->col(group_B->at(B)) );
                    if (dist_new < lower_cutoff) {
                        return PAIR_INFINITY;
                    }
                    else if (dist_new < upper_cutoff) {
                        Delta_E += Potential(dist_new);
                    }

                    dist_old = arma::norm( pos_old->col(group_A->at(Ap)) - pos_old->col(group_B->at(B)) );
                    if (dist_old < upper_cutoff) {
                        Delta_E -= Potential(dist_old);
                    }
                }
            }
        }

        /*
            Moved - Moved checks
        */
        for (int mov2=mov+1;mov2<num_moved;mov2++){
            ref1 = mov2*2;
            ref2 = ref1+1;
            A1p = A_moved->at(ref1);
            A2p = A_moved->at(ref2);
            B1p = B_moved->at(ref1);
            B2p = B_moved->at(ref2);

            // check A with Bp
            for (int A=A1;A<=A2;A++) {
                for (int Bp=B1p;Bp<=B2p;Bp++) {
                    dist_new = arma::norm( pos_new->col(group_A->at(A)) - pos_new->col(group_B->at(Bp)) );
                    if (dist_new < lower_cutoff) {
                        return PAIR_INFINITY;
                    }
                    else if (dist_new < upper_cutoff) {
                        Delta_E += Potential(dist_new);
                    }

                    dist_old = arma::norm( pos_old->col(group_A->at(A)) - pos_old->col(group_B->at(Bp)) );
                    if (dist_old < upper_cutoff) {
                        Delta_E -= Potential(dist_old);
                    }
                }
            }

            // check B with Ap
            for (int B=B1;B<=B2;B++) {
                for (int Ap=A1p;Ap<=A2p;Ap++) {
                    dist_new = arma::norm( pos_new->col(group_A->at(Ap)) - pos_new->col(group_B->at(B)) );
                    if (dist_new < lower_cutoff) {
                        return PAIR_INFINITY;
                    }
                    else if (dist_new < upper_cutoff) {
                        Delta_E += Potential(dist_new);
                    }

                    dist_old = arma::norm( pos_old->col(group_A->at(Ap)) - pos_old->col(group_B->at(B)) );
                    if (dist_old < upper_cutoff) {
                        Delta_E -= Potential(dist_old);
                    }
                }
            }
        }
    }
    /*
        Individual checks
    */

    if (single_group) {
        for (int ind=0;ind<num_indiv;ind++) {
            ref1 = ind*2;
            ref2 = ref1+1;
            A1 = A_individual->at(ref1);
            A2 = A_individual->at(ref2);

            for (int A=A1;A<A2;A++) {
                for (int Ap=A+1;Ap<=A2;Ap++) {

                    dist_new = arma::norm( pos_new->col(group_A->at(A)) - pos_new->col(group_A->at(Ap)) );
                    if (dist_new < lower_cutoff) {
                        return PAIR_INFINITY;
                    }
                    else if (dist_new < upper_cutoff) {
                        Delta_E += Potential(dist_new);
                    }

                    dist_old = arma::norm( pos_old->col(group_A->at(A)) - pos_old->col(group_A->at(Ap)) );
                    if (dist_old < upper_cutoff) {
                        Delta_E -= Potential(dist_old);
                    }
                }
            }
        }
    }
    else {
        for (int ind=0;ind<num_indiv;ind++) {
            ref1 = ind*2;
            ref2 = ref1+1;
            A1 = A_individual->at(ref1);
            A2 = A_individual->at(ref2);
            B1 = B_individual->at(ref1);
            B2 = B_individual->at(ref2);

            // check A with B
            for (int A=A1;A<=A2;A++) {
                for (int B=B1;B<=B2;B++) {
                    dist_new = arma::norm( pos_new->col(group_A->at(A)) - pos_new->col(group_B->at(B)) );
                    if (dist_new < lower_cutoff) {
                        return PAIR_INFINITY;
                    }
                    else if (dist_new < upper_cutoff) {
                        Delta_E += Potential(dist_new);
                    }

                    dist_old = arma::norm( pos_old->col(group_A->at(A)) - pos_old->col(group_B->at(B)) );
                    if (dist_old < upper_cutoff) {
                        Delta_E -= Potential(dist_old);
                    }
                }
            }
        }
    }
    return Delta_E;
}


/*---------------------------------------------------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------------------------------------------------*/


/* ------------------------------------------------------------------------- */


double Pair::Potential(double distance) {
    std::cout << "Parental" << std::endl;
    std::exit(0);
    return 1e10;
}

unsigned Pair::find_id_in_group_A(int monomer_id) {
    for (unsigned i=0;i<group_A->size();i++) {
        if (group_A->at(i) == monomer_id) {
            return i;
        }
    }
    return -1;
}

unsigned Pair::find_id_in_group_B(int monomer_id) {
    for (unsigned i=0;i<group_B->size();i++) {
        if (group_B->at(i) == monomer_id) {
            return i;
        }
    }
    return -1;
}






