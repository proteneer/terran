#include "MethodsGaussian.h"
#include <stdexcept>
#include <assert.h>
#include <iostream>
#include <fstream>
#include <algorithm>  

using namespace std;

namespace Terran {

// useful utility for plotting the gaussin
static void plotGaussian(const vector<Param> &params) {
    ofstream mixture("mixture.dat");
    for(double xn = -10; xn < 10; xn += 0.01) {
        mixture << xn << " " << gaussianMixture(params, xn) << endl;
    }
    ofstream mixtureDx("mixtureDx.dat");
    for(double xn = -10; xn < 10; xn += 0.01) {
        mixtureDx << xn << " " << gaussianMixtureDx(params, xn) << endl;
    }
    ofstream mixtureDx2("mixtureDx2.dat");
    for(double xn = -10; xn < 10; xn += 0.01) {
        mixtureDx2 << xn << " " << gaussianMixtureDx2(params, xn) << endl;
    }
}

static bool compare_by_x(double2 a, double2 b) {
	return a.x < b.x;
}

MethodsGaussian::MethodsGaussian(const vector<Param> &params) : Methods(params) {
	
	// compute left and right bounds
    double min =  numeric_limits<double>::max();
    double max = -numeric_limits<double>::max();
	double max_sig = -numeric_limits<double>::max();
    for(int i=0; i < params_.size(); i++) {
        min = (params_[i].u < min) ? params_[i].u : min;
        max = (params_[i].u > max) ? params_[i].u : max;
        max_sig = (params_[i].s > max) ? params_[i].s : max;
    }
 
    while(gaussianMixture(params_,min) > 1e-5)
        min -= max_sig;
    while(gaussianMixture(params_,max) > 1e-5)
        max += max_sig;
 
	// sample from the distribution directly
	vector<double2> samples;
	for(int i=0; i <2500; i++) {
		double2 sample;
		sample.x = gaussianMixtureSample(params_);
		sample.y = gaussianMixture(params_, sample.x);
		samples.push_back(sample);
	}

	sort(samples.begin(), samples.end(), compare_by_x);

	// scan through the values and identify brackets
	for(int i=0; i < samples.size(); i++) {
		if(i !=0 && i != samples.size()-1) {
			double2 left = samples[i-1];
			double2 middle = samples[i];
			double2 right = samples[i+1];
		
			if(middle.y > left.y && middle.y > right.y) {
				maxBrackets_.push_back(Bracket(left.x, middle.x, right.x));
			} else if(middle.y < left.y && middle.y < right.y) {
				minBrackets_.push_back(Bracket(left.x, middle.x, right.x));
			}
		}
	}

	if(maxBrackets_.size() == 0) {
		throw(std::runtime_error("MethodsGaussian::MethodsGaussian() - maxBrackets_ is empty"));
	}

}

// golden section search, from numerical recipes
static inline void shift3(double &a, double &b, double &c, double d) {
	a=b;
	b=c;
	c=d;
}

static inline void shift2(double &a, double &b, double c) {
	a=b;
	b=c;
}

// This implementation uses the golden section algorithm from numerical recipes.
// The initial brackets are computed by sampling points from the initial distribution.
vector<double> MethodsGaussian::findMaxima() const {
	vector<double> maxima;

	const double R = 0.61803399;
	const double C = 1.0 - R;
	const double tolerance = 1e-7;

	for(int i=0; i < maxBrackets_.size(); i++) {	
		double ax = maxBrackets_[i].left;
		double bx = maxBrackets_[i].middle;
		double cx = maxBrackets_[i].right;
		double f1,f2,x0,x1,x2,x3;
		x0 = ax;
		x3 = cx;

		if(fabs(cx-bx) > fabs(bx-ax)) {
			x1 = bx;
			x2 = bx + C*(cx-bx);
		} else{ 
			x2 = bx;
			x1 = bx - C*(bx-ax);
		}
		f1 = gaussianMixture(params_, x1);
		f2 = gaussianMixture(params_, x2);

		while(fabs(x3-x0) > tolerance*(fabs(x1)+fabs(x2))) {
			if(f2 < f1) {
				shift3(x3,x2,x1,R*x2+C*x0);
				shift2(f2,f1,gaussianMixture(params_,x1));
			} else {
				shift3(x0,x1,x2,R*x1+C*x3);
				shift2(f1,f2,gaussianMixture(params_,x2));
			}
		}
		maxima.push_back(x1);
	}
	return maxima;
}


// TODO: This does redundant work, can adapt findMaxima() to do this.
vector<double> MethodsGaussian::findMinima() const {
    vector<double> minima;
    vector<double> maxima = findMaxima();
    assert(maxima.size() > 0);
    sort(maxima.begin(), maxima.end());
    const double perturb = 5e-2;

    if(maxima.size() > 1) {
        for(int k=0; k<maxima.size()-1; k++) {
            double ak = maxima[k];
            double bk = maxima[k+1];
            ak += perturb;
            bk -= perturb;
            double mk = 0;
            double my = 0;
            int iteration = 0;
            // bisect the interval
            do {
                if(iteration > 1e5) {
                    throw(std::runtime_error("MethodsGaussian::findMinima() - maximized number of iterations reached."));
                }
                else
                    iteration++;
                
                mk = (ak+bk)/2.0;
                my = gaussianMixtureDx(params_, mk);
                if(my < 0) 
                    ak = mk;
                else
                    bk = mk;
            } while (fabs(ak-bk) > 1e-8);
            minima.push_back(mk);
        }
    }
    return minima;
}

} // namespace Terran
