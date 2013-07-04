#ifndef MATHFUNCTIONS_H
#define MATHFUNCTIONS_H

#include <math.h>
#include <vector>
#include "Param.h"

const double PI = 3.14159265358;

inline double gaussian(double uk, double sk, double xn) {
    return 1.0/(sqrt(2*PI)*sk)*exp(-(0.5)*pow((xn-uk)/sk,2));
} 

inline double periodicGaussian(double uk, double sk, double xn, int numImages, double period) {
    double sum = 0;
    for(int r=-numImages; r<=numImages; r++) {
        sum += gaussian(uk+r*period, sk, xn);
    }
    return sum;
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

inline double periodicGaussianMixture(const std::vector<Param> &params, double xn, int numImages = 20, int period = 2*PI) {
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

// TODO: Derivatives of both GMs and PGMs
#endif
