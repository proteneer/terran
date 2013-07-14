#ifndef MATHFUNCTIONS_H
#define MATHFUNCTIONS_H

#include <math.h>
#include <vector>
#include "Param.h"

const double PI = 3.14159265358;

inline double gaussian(double uk, double sk, double xn) {
    return 1.0/(sqrt(2*PI)*sk)*exp(-(0.5)*pow((xn-uk)/sk,2));
} 

inline double gaussianDx(double uk, double sk, double xn) {
    double multiplier = (uk-xn)/(sk*sk);
    double suffix = multiplier*gaussian(uk, sk, xn);
    return multiplier*suffix;
}

inline double periodicGaussian(double uk, double sk, double xn, int numImages, double period) {
    double sum = 0;
    for(int r=-numImages; r<=numImages; r++) {
    // this should be gaussian(uk-r*period) to match the paper since r is symmetric.
        sum += gaussian(uk-r*period, sk, xn);
    }
    return sum;
}


/*
Equation [8] in the paper is completely wrong!
*/

inline double periodicGaussianDx(double uk, double sk, double xn, int numImages, double period) {
    double prefactor = 1/(sqrt(2*PI)*sk*sk*sk);
    double sum = 0;
    for(int r=-numImages; r<=numImages; r++) {
        sum += (uk-xn+r*period)*exp(-0.5*pow((xn-uk-r*period)/sk,2));
    }
    return prefactor*sum;
}

inline double gaussianMixture(const std::vector<Param> &params, double xn) {
    double sum = 0;
    for(int k=0; k<params.size(); k++) {
        double pk = params[k].p;
        double uk = params[k].u;
        double sk = params[k].s;
        sum += pk*gaussian(uk, sk, xn);
    }
    return sum;
}

inline double gaussianMixtureDx(const std::vector<Param> &params, double xn) {
    double sum = 0;
    for(int k=0; k<params.size(); k++) {
        double pk = params[k].p;
        double uk = params[k].u;
        double sk = params[k].s;
        sum += pk*gaussianDx(uk, sk, xn);
    }
    return sum;

}

inline double periodicGaussianMixture(const std::vector<Param> &params, double xn, double period = 2*PI, int numImages = 10) {
    //assert(xn >= (-period/2-1e-6) && xn <= (period/2+1e-6));  
    double sum = 0;
    for(int k=0; k<params.size(); k++) {
        double pk = params[k].p;
        double uk = params[k].u;
        double sk = params[k].s;
        sum += pk*periodicGaussian(uk, sk, xn, numImages, period);
    }
    return sum;
}

inline double periodicGaussianMixtureDx(const std::vector<Param> &params, double xn, double period = 2*PI, int numImages = 10) {
    //assert(xn >= (-period/2-1e-6) && xn <= (period/2+1e-6));  
    double sum = 0;
    for(int k=0; k<params.size(); k++) {
        double pk = params[k].p;
        double uk = params[k].u;
        double sk = params[k].s;
        sum += pk*periodicGaussianDx(uk, sk, xn, numImages, period);
    }
    return sum;
}

// TODO: Derivatives of both GMs and PGMs
#endif
