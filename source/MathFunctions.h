#ifndef MATHFUNCTIONS_H
#define MATHFUNCTIONS_H

#include <math.h>
#include <vector>
#include "Param.h"


#include <iostream>
using namespace std;

const double PI = 3.14159265358;

namespace Terran {

// moves a point to inside the periodic domain
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

// evaluates the gaussian function at point xn given parameters uk and sk
inline double gaussian(double uk, double sk, double xn) {
    return 1.0/(sqrt(2*PI)*sk)*exp(-(0.5)*pow((xn-uk)/sk,2));
} 

// derivative of the gaussian
inline double gaussianDx(double uk, double sk, double xn) {
    double multiplier = (uk-xn)/(sk*sk);
    double suffix = multiplier*gaussian(uk, sk, xn);
    return multiplier*suffix;
}

// evaluates the periodic gaussian function at point xn
// uk and sk retain the semantic of mean and standard deviation
// numImages is the number of times to wrap around
// period is defined by the domain of the function

// this version of the periodicGaussian picks the number of images automatically based 
// the value of sk and uk

/*
inline double periodicGaussian(double uk, double sk, double xn, double period) {
    double sum = 0;


    for(int r=-numImages; r<=numImages; r++) {
        sum += gaussian(uk-r*period, sk, xn);
    }
    return sum;
}
*/

inline double periodicGaussian(double uk, double sk, double xn, int numImages, double period) {
    double sum = 0;
    for(int r=-numImages; r<=numImages; r++) {
        sum += gaussian(uk-r*period, sk, xn);
    }
    return sum;
}

/*
Equation [8] in the paper is completely wrong!
*/



// derivative of the periodic gaussian function


//   |----0----|
// --|u---0----|--
// l-|u-r-0----|--


//     **
//     **
//   ******
// **********


//                     |               l  |  u     r         |
// -----------------------------------------------------------------------------
// -2-2-2-2-2-2-1-1-1-1-1-1-1-1-1-1000000000000000001111111111111111111122222222              p

//                        
//                     |               l  |  u     r         |
// -----------------------------------------------------------------------------
// -2-2-2-2-2-2-2-2-2-2-1-1-1-1-1-1-1-1-1 00000000000000000001111111111111111111

inline double periodicGaussian(double uk, double sk, double xn, double period) {
    // six standard deviations covers 99.99966% of the total probability
    double left = uk-7*sk;
    double right = uk+7*sk;


    int lower = floor(left/period+0.5);
    int upper = floor(right/period+0.5);
    
    // left/period 

    //lower = -2;
    //upper = 2;

    cout << left << " " << right << " " << lower << " " << upper << endl;

    int numImages = 10;
    double sum = 0;
    for(int r=lower; r<=upper; r++) {
        sum += gaussian(uk-r*period, sk, xn);
    }
    return sum;
}

// derivative of the periodic gaussian function
inline double periodicGaussianDx(double uk, double sk, double xn, int numImages, double period) {
    double prefactor = 1/(sqrt(2*PI)*sk*sk*sk);
    double sum = 0;
    for(int r=-numImages; r<=numImages; r++) {
        sum += (uk-xn+r*period)*exp(-0.5*pow((xn-uk-r*period)/sk,2));
    }
    return prefactor*sum;
}

// evaluates xn given a gaussian mixture
// where sum pk = 1
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

// derivative of a gaussian mixture
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

// evaluates xn given a periodic gaussian mixture
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

// derivative of a gaussian mixture
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
        double probability = periodicGaussian(u, s, x1, 100, period);
        // make a random sample between 0 and 1 
        double draw = (double) rand() / (double) RAND_MAX;
        if(draw < probability) {
            return x1;
        }
    }
}

}

#endif
