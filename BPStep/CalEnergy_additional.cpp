#include "BPStep.h"
#ifndef BPS_USE_EVALENERGY

bool BPStep::cal_MC_isotropic_bending() {

    std::cout << "isotropic_bending starting" << std::endl;

    std::random_device                      rd{};
    std::mt19937                            gen{rd()};
    std::uniform_real_distribution<double>  uniform{0.0,2*M_PI};
    std::seed_seq seed{rd(), rd(), rd(), rd(), rd(), rd(), rd(), rd()};
    gen.seed(seed);

    double E1,E2;
    double phi,cosphi,sinphi;
    arma::colvec Theta1,Theta2;
    Theta2 = arma::zeros(3);
    for (long long int i=0;i<BPS_ISOTROPIC_BENDING_NUM_SAMPLED;i++){
        Theta1    = {uniform(gen),uniform(gen),uniform(gen)};
        phi      = uniform(gen);
        sinphi   = std::sin(phi);
        cosphi   = std::cos(phi);
        E1       = cal_beta_energy_diag(Theta1);
        Theta2(0) = Theta1(0)*cosphi-Theta1(1)*sinphi;
        Theta2(1) = Theta1(0)*sinphi+Theta1(1)*cosphi;
        Theta2(2) = Theta1(2);
        E2       = cal_beta_energy_diag(Theta2);
        if (std::abs(E2-E1)>BPS_ISOTROPIC_BENDING_EQUAL_THRES) {
            std::cout << std::abs(E2-E1) << std::endl;
            std::cout << "Is NOT isotropic!" << std::endl;
            return false;
        }
    }
    std::cout << "Is isotropic!" << std::endl;
    return true;
}

////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

arma::mat BPStep::cal_MC_cov(long long int steps,double sigma) {
    std::random_device                      rd{};
    std::mt19937                            gen{rd()};
    std::normal_distribution<double>        normaldist{0.0,1.0};
    std::uniform_real_distribution<double>  uniformdist{0.0,1.0};
    std::seed_seq seed{rd(), rd(), rd(), rd(), rd(), rd(), rd(), rd()};
    gen.seed(seed);

    arma::mat    cov    = arma::zeros(3,3);
    arma::colvec Theta    = arma::zeros(3);
    arma::colvec newTheta = arma::zeros(3);
    double currentbetaE,newbetaE;
    currentbetaE = cal_beta_energy_diag(Theta);
    for (long long int i=0;i<steps;i++){
        newTheta(0) = Theta(0) + normaldist(gen)*sigma;
        newTheta(1) = Theta(1) + normaldist(gen)*sigma;
        newTheta(2) = Theta(2) + normaldist(gen)*sigma;
        newbetaE = cal_beta_energy_diag(newTheta);
        if (std::exp( currentbetaE - newbetaE ) > uniformdist(gen)) {
            Theta         = newTheta;
            currentbetaE  = newbetaE;
        }
        cov    += Theta*Theta.t();
    }
    return cov /= steps;
}

#endif
