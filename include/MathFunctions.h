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

struct double2 {
	double x;
	double y;
};

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

// ------------------
// Canonical Gaussian
// ------------------
//
// The following 4 functions expression define gaussian, gaussian mixtures, and their derivatives

inline double gaussian(double uk, double sk, double xn) {
    return 1.0/(sqrt(2*PI)*sk)*exp(-(0.5)*pow((xn-uk)/sk,2));
} 

// These derivatives are not multiplied sqrt(2PI)
inline double gaussianDx(double uk, double sk, double xn) {
    double multiplier = (uk-xn)/(sk*sk);
    double suffix = gaussian(uk, sk, xn);
    return multiplier*suffix;
}

inline double gaussianDx2(double uk, double sk, double xn) {
    double multiplier = ((uk-xn)*(uk-xn)-sk*sk)/(sk*sk*sk*sk);
    double suffix = gaussian(uk, sk, xn);
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

inline double gaussianMixtureDx2(const std::vector<Param> &params, double xn) {
    double sum = 0;
    for(int k=0; k<params.size(); k++) {
        double pk = params[k].p;
        double uk = params[k].u;
        double sk = params[k].s;
        sum += pk*gaussianDx2(uk, sk, xn);
    }
    return sum;
}

// -------------------------------------------
// Periodic Gaussian w/ automatic image tuning
// -------------------------------------------
//
// Periodic gaussians are basically wrapped versions of the canonical
// gaussian. These functions below guess the num of images needed to
// approximate the gaussian to a very high degree of accuracy (1e-9)

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
    double multiplier = 1/(sk*sk);
    // heuristic derived from tuning
    int numImages = ceil(sk)+1;
    double sum = 0;
    for(int r=-numImages; r<=numImages; r++) {
        sum += (uk-xn+r*period)*gaussian(uk+r*period, sk, xn);
    }
    return multiplier*sum;
}

inline double periodicGaussianDx2(double uk, double sk, double xn, double period) {
    double multiplier = 1/(sk*sk*sk*sk);
    double sum = 0;
    int numImages = 10;
    for(int r=-numImages; r<=numImages; r++) {
        sum += ((uk-xn+r*period)*(uk-xn+r*period)-sk*sk)*gaussian(uk+r*period, sk, xn);
    }
    return sum*multiplier;
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

inline double periodicGaussianMixtureDx2(const std::vector<Param> &params, double xn, double period) {
    //assert(xn >= (-period/2-1e-6) && xn <= (period/2+1e-6));  
    double sum = 0;
    for(int k=0; k<params.size(); k++) {
        double pk = params[k].p;
        double uk = params[k].u;
        double sk = params[k].s;
        sum += pk*periodicGaussianDx2(uk, sk, xn, period);
    }
    return sum;
}

// ----------------------------------------
// Periodic Gaussian w/ manual image tuning
// ----------------------------------------
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

// ------------------------------------------
// Methods for sampling from 1D distributions
// ------------------------------------------

// draws a random sample from the periodic gaussian with mean u,
// std. deviation s, and period period
// assume x is in -period/2, period/2
inline double periodicGaussianSample(double u, double s, double period) {

	if(s < 0.0002) {
		throw(std::runtime_error("s is too small"));
	}

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
			//cout << "good" << endl;
            return x1;
        } else {
		
			//cout << u << " " << s << period << " " << draw << " " << probability <<  endl;
		
		}
    }
}

// draws a random sample from the periodic gaussian with mean u,
// std. deviation s
inline double gaussianSample(double u, double s) {
    
    double z1 = (double)rand()/(double)RAND_MAX;
    double z2 = (double)rand()/(double)RAND_MAX;

    double x = sqrt(-2*log(z1))*cos(2*PI*z2);
    x *= s;
    x += u;
    return x;
}

// draw a random sample from the mixture params
inline double gaussianMixtureSample(vector<Param> params) {

	// running sum of component probabilities
	// (this can be cached if need be)
	double cumulant = 0;
	vector<double> running_sum(params.size());
	for(int i=0; i < params.size(); i++) {
		cumulant += params[i].p;
		running_sum[i] = cumulant;
	}
	
	double p1 = (double) rand() / (double) RAND_MAX;

	int k = 0;
	while(running_sum[k] < p1) {
		k++;
	}

	return gaussianSample(params[k].u, params[k].s);
}

inline double periodicGaussianMixtureSample(vector<Param> params, double period) {

	// running sum of component probabilities
	// (this can be cached if need be)
	double cumulant = 0;
	vector<double> running_sum(params.size());
	for(int i=0; i < params.size(); i++) {
		cumulant += params[i].p;
		running_sum[i] = cumulant;
	}
	
	double p1 = (double) rand() / (double) RAND_MAX;

	if(p1 > 1) {
		cout << p1 << endl;
		throw(std::runtime_error("p1 > 1"));
	}

	int k = 0;
	
	for(; k < running_sum.size(); k++) 
		if(p1 <= running_sum[k]) {
			break;
		}

	if(k >= running_sum.size()) {
		k = running_sum.size()-1;
	}
	//cout << k << ":" << params[k].u << " " << params[k].s << endl;

	if(params[k].s < 0.001) {
		cout << p1 << endl;
		cout << params.size() << endl;
		for(int i=0; i < running_sum.size(); i++) {
			cout << running_sum[i] << endl;
		}
		cout << "-------------------" << endl;
		int a = 0;
		for(; a < running_sum.size(); a++) {
			cout << p1 << " " << a << " " << running_sum[a] << endl;
			if(p1 <= running_sum[a]) {
				cout << "breaking a" << endl;
				break;
			}
		}

		if(p1 > running_sum[6]) {
			cout.precision(20);
			cout << std::fixed << p1 << endl;
			cout << std::fixed << running_sum[6] << endl;
			cout << "WTF" << endl;
		}

		cout << "A: " << a << endl;

	}

	return periodicGaussianSample(params[k].u, params[k].s, period);
}

}// namespace Terran

#endif
