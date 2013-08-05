#ifndef MATHFUNCTIONS_H
#define MATHFUNCTIONS_H

#include <math.h>
#include <vector>
#include "Param.h"

#ifdef _WINDOWS
#define isnan(x) _isnan(x) 
#endif

#include <iostream>
using namespace std;

const double PI = 3.14159265358;

namespace Terran {

// moves a point to inside the periodic domain
// default arguments are dangerous!
inline double normalize(double x, double left=-PI, double right=PI) {
    double period = right-left;
    while(x > right) {
        x -= period;
    }
    while(x < left) {
        x += period;
    }
    return x;
}

// calculates the shortest distance between two points given a period
inline double periodicDifference(double x1, double x2, double period) {
    double diff = x1-x2;
    diff -= floor(diff/period+0.5)*period;
    return diff;
}

// returns |x1-x2| given shortest period
inline double fabsp(double x1, double x2, double period) {
    return fabs(periodicDifference(x1,x2,period));
}

// ******************
// Canonical Gaussian
// ******************
//
// The following 4 functions expression define gaussian, gaussian mixtures, and their derivatives

inline double gaussian(double uk, double sk, double xn) {
    return 1.0/(sqrt(2*PI)*sk)*exp(-(0.5)*pow((xn-uk)/sk,2));
} 

inline double gaussianDx(double uk, double sk, double xn) {
    double multiplier = (uk-xn)/(sk*sk);
    double suffix = multiplier*gaussian(uk, sk, xn);
    return multiplier*suffix;
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

// *******************************************
// Periodic Gaussian w/ automatic image tuning
// *******************************************
//
// Periodic gaussians are basically wrapped versions of the canonical
// gaussian. These functions below guess the num of images needed to
// approximate the gaussian to a very high degree of accuracy (1e-9)

/*
Equation [8] in the paper is completely wrong!
*/

inline double periodicGaussian(double uk, double sk, double xn, double period) {
    double left = uk-7*sk;
    double right = uk+7*sk;
    int lower = floor(left/period+0.5);
    int upper = floor(right/period+0.5);
    double sum = 0;
    for(int r=lower; r<=upper; r++) {
        sum += gaussian(uk-r*period, sk, xn);
    }
    return sum;
}


inline double periodicGaussianDx(double uk, double sk, double xn, double period) {
    double prefactor = 1/(sqrt(2*PI)*sk*sk*sk);
    // heuristic derived from tuning
    int numImages = ceil(sk)+1;
    double sum = 0;
    for(int r=-numImages; r<=numImages; r++) {
        sum += (uk-xn+r*period)*exp(-0.5*pow((xn-uk-r*period)/sk,2));
    }
    return prefactor*sum;
}

inline double periodicGaussianMixture(const std::vector<Param> &params, double xn, double period) {
    //assert(xn >= (-period/2-1e-6) && xn <= (period/2+1e-6));  
    double sum = 0;
    for(int k=0; k<params.size(); k++) {
        double pk = params[k].p;
        double uk = params[k].u;
        double sk = params[k].s;
        sum += pk*periodicGaussian(uk, sk, xn, period);
    }
    return sum;
}

inline double periodicGaussianMixtureDx(const std::vector<Param> &params, double xn, double period) {
    //assert(xn >= (-period/2-1e-6) && xn <= (period/2+1e-6));  
    double sum = 0;
    for(int k=0; k<params.size(); k++) {
        double pk = params[k].p;
        double uk = params[k].u;
        double sk = params[k].s;
        sum += pk*periodicGaussianDx(uk, sk, xn, period);
    }
    return sum;
}

// ****************************************
// Periodic Gaussian w/ manual image tuning
// ****************************************
// 
// These overloaded versions of the periodic gaussian allow for finer control over the number of images
// used. These functions are used in unit tests

inline double periodicGaussian(double uk, double sk, double xn, double period, int numImages) {
    double sum = 0;
    for(int r=-numImages; r<=numImages; r++) {
        sum += gaussian(uk-r*period, sk, xn);
    }
    return sum;
}

inline double periodicGaussianDx(double uk, double sk, double xn, double period, int numImages) {
    double prefactor = 1/(sqrt(2*PI)*sk*sk*sk);
    double sum = 0;
    for(int r=-numImages; r<=numImages; r++) {
        sum += (uk-xn+r*period)*exp(-0.5*pow((xn-uk-r*period)/sk,2));
    }
    return prefactor*sum;
}

inline double periodicGaussianMixture(const std::vector<Param> &params, double xn, double period, int numImages) {
    //assert(xn >= (-period/2-1e-6) && xn <= (period/2+1e-6));  
    double sum = 0;
    for(int k=0; k<params.size(); k++) {
        double pk = params[k].p;
        double uk = params[k].u;
        double sk = params[k].s;
        sum += pk*periodicGaussian(uk, sk, xn, period, numImages);
    }
    return sum;
}

inline double periodicGaussianMixtureDx(const std::vector<Param> &params, double xn, double period, int numImages) {
    //assert(xn >= (-period/2-1e-6) && xn <= (period/2+1e-6));  
    double sum = 0;
    for(int k=0; k<params.size(); k++) {
        double pk = params[k].p;
        double uk = params[k].u;
        double sk = params[k].s;
        sum += pk*periodicGaussianDx(uk, sk, xn, period, numImages);
    }
    return sum;
}

// draws a random sample from the periodic gaussian with mean u,
// std. deviation s, and period period
// assume x is in -period/2, period/2
inline double periodicGaussianSample(double u, double s, double period) {
    // sample between [0,1]
    while(true) {
        double x1 = (double) rand() / (double) RAND_MAX;
        // move to [-0.5, 0.5]
        x1 -= 0.5;
        // move to [-0.5*2PI, -0.5*2PI]
        x1 *= period;
        // determine the probability of drawing x1
        double probability = periodicGaussian(u, s, x1, period);
        // make a random sample between 0 and 1 
        double draw = (double) rand() / (double) RAND_MAX;
        if(draw < probability) {
            return x1;
        }
    }
}

}

#endif
