#include "ExtraArma.h"

int arma_argmax(arma::colvec vec)
{
    double max_val = vec(0);
    int    max_id  = 0;
    for (int i=1;i<vec.n_elem;i++) {
        if (vec(i) > max_val) {
            max_val = vec(i);
            max_id  = i;
        }
    }
    return max_id;
}

arma::colvec normalize(arma::colvec vec) {
    double norm = arma::norm(vec);
    if (std::abs(norm) < 1e-14) {
        return arma::zeros(3);
    }
    return vec/norm;
}
