#ifndef __SO3METHODS_INCLUDED__
#define __SO3METHODS_INCLUDED__

#include <armadillo> // armadillo
#include <iostream>  // input output
#include <stdexcept> // exceptions
#include <cmath>     //
#include <string>    // string
#include <algorithm>


#ifndef __SO3METHODS_CONSTANTS__
#define __SO3METHODS_CONSTANTS__

#define SO3M_EPSILON                1e-12
#define S03M_CLOSE_TO_ONE           0.999999999999
#define S03M_CLOSE_TO_MINUS_ONE    -0.999999999999

#endif

//arma::mat       getRotMat(const arma::colvec& Omega);
//arma::colvec    ExtractTheta(const arma::mat& R);
//double          ExtractTheta1(const arma::mat& R);
//double          ExtractTheta2(const arma::mat& R);
//double          ExtractTheta3(const arma::mat& R);
//arma::mat       Rotz(const double& theta);

inline arma::mat getRotMat(const arma::colvec& Omega)
{
    double Om = arma::norm(Omega);
    if (Om < SO3M_EPSILON) {
        return arma::eye(3,3);
    }
    double cosOm = std::cos(Om);
    double sinOm = std::sin(Om);
    double Omsq	 = Om*Om;
    double fac1  = (1-cosOm)/Omsq;
    double fac2  = sinOm/Om;
    arma::mat R  = arma::zeros(3,3);
    R(0,0) = cosOm+Omega(0)*Omega(0)*fac1;
    R(1,1) = cosOm+Omega(1)*Omega(1)*fac1;
    R(2,2) = cosOm+Omega(2)*Omega(2)*fac1;
    double A = Omega(0)*Omega(1)*fac1;
    double B = Omega(2)*fac2;
    R(0,1) = A-B;
    R(1,0) = A+B;
    A = Omega(0)*Omega(2)*fac1;
    B = Omega(1)*fac2;
    R(0,2) = A+B;
    R(2,0) = A-B;
    A = Omega(1)*Omega(2)*fac1;
    B = Omega(0)*fac2;
    R(1,2) = A-B;
    R(2,1) = A+B;
    return R;
}

inline arma::colvec ExtractTheta(const arma::mat& R) {
    double val = 0.5* (arma::trace(R)-1);
    if (val > S03M_CLOSE_TO_ONE) {
        return arma::zeros(3);
    }
    if (val < S03M_CLOSE_TO_MINUS_ONE) {
        if (R(0,0) > S03M_CLOSE_TO_ONE) {
            return {M_PI,0,0};
        }
        if (R(1,1) > S03M_CLOSE_TO_ONE) {
            return {0,M_PI,0};
        }
        return {0,0,M_PI};
    }
    double Th = std::acos(val);
    arma::colvec Theta =  {(R(2,1)-R(1,2)),(R(0,2)-R(2,0)),(R(1,0)-R(0,1))};
    Theta = Th*0.5/std::sin(Th) * Theta;
    return Theta;
}

inline double ExtractTheta1(const arma::mat& R) {
    double val = 0.5* (arma::trace(R)-1);
    if (val > S03M_CLOSE_TO_ONE) {
        return 0;
    }
    if (val < S03M_CLOSE_TO_MINUS_ONE) {
        if (R(0,0) > S03M_CLOSE_TO_ONE) {
            return M_PI;
        }
        return 0;
    }
    double Th = std::acos(val);
    return Th*0.5/std::sin(Th) * (R(2,1)-R(1,2));
}

inline double ExtractTheta2(const arma::mat& R) {
    double val = 0.5* (arma::trace(R)-1);
    if (val > S03M_CLOSE_TO_ONE) {
        return 0;
    }
    if (val < S03M_CLOSE_TO_MINUS_ONE) {
        if (R(1,1) > S03M_CLOSE_TO_ONE) {
            return M_PI;
        }
        return 0;
    }
    double Th = std::acos(val);
    return Th*0.5/std::sin(Th) * (R(0,2)-R(2,0));
}

inline double ExtractTheta3(const arma::mat& R) {
    double val = 0.5* (arma::trace(R)-1);
    if (val > S03M_CLOSE_TO_ONE) {
        return 0;
    }
    if (val < S03M_CLOSE_TO_MINUS_ONE) {
        if (R(2,2) > S03M_CLOSE_TO_ONE) {
            return M_PI;
        }
        return 0;
    }
    double Th = std::acos(val);
    return Th*0.5/std::sin(Th)* (R(1,0)-R(0,1));
}

inline arma::mat Rotz(const double& theta) {
    double c = std::cos(theta);
    double s = std::sin(theta);
    arma::mat Rz = {{c,-s,0},{s,c,0},{0,0,1}};
    return Rz;
}

#endif
