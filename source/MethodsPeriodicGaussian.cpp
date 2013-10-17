#include "MethodsPeriodicGaussian.h"
#include <stdexcept>
#include <assert.h>
#include <iostream>
#include <fstream>
#include <algorithm>  
#include <sstream>

using namespace std;


#include "omp.h"


namespace Terran {

static void plotPeriodicGaussian(const vector<Param> &params, double period, int images) {
    ofstream mixture("mixture.dat");
    for(double xn = -period/2; xn < period/2; xn += 0.01) {
        mixture << xn << " " << periodicGaussianMixture(params, xn, period, images) << endl;
    }
    ofstream mixtureDx("mixtureDx.dat");
    for(double xn = -period/2; xn < period/2; xn += 0.01) {
        mixtureDx << xn << " " << periodicGaussianMixtureDx(params, xn, period, images) << endl;
    }
}

static bool compare_by_x(double2 a, double2 b) {
	return a.x < b.x;
}

MethodsPeriodicGaussian::MethodsPeriodicGaussian(const vector<Param> &params, 
    double period) : Methods(params), period_(period) {
	// sample from the distribution directly
	vector<double2> samples;
	for(int i=0; i <2500; i++) {
		double2 sample;
		sample.x = periodicGaussianMixtureSample(params_, period_);
		sample.y = periodicGaussianMixture(params_, sample.x, period_);
		samples.push_back(sample);
	}
	sort(samples.begin(), samples.end(), compare_by_x);
	vector<double2> samples2;
	for(int i=1; i < samples.size(); i++) {
		if(samples[i].x != samples[i-1].x) {
			samples2.push_back(samples[i]);
		}
	}
	samples = samples2;
	// scan through the values and identify potential brackets cyclically
	for(int i=0; i < samples.size(); i++) {
		double2 left = samples[(i-1+samples.size())%samples.size()];
		double2 middle = samples[i];
		double2 right = samples[(i+1+samples.size())%samples.size()];
		if((middle.y > left.y) && (middle.y > right.y)) {
			maxBrackets_.push_back(Bracket(left.x, middle.x, right.x));
		}
	}

	if(maxBrackets_.size() == 0) {
		throw(std::runtime_error("MethodsPeriodicGaussian::MethodsPeriodicGaussian() - maxBrackets_ is empty"));
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

vector<double> MethodsPeriodicGaussian::findMaxima() const {
	vector<double> maxima;
	const double R = 0.61803399;
	const double C = 1.0 - R;
	const double tolerance = 1e-7;

	for(int i=0; i < maxBrackets_.size(); i++) {
		double ax = maxBrackets_[i].left;
		double bx = maxBrackets_[i].middle;
		double cx = maxBrackets_[i].right;
	}

	for(int i=0; i < maxBrackets_.size(); i++) {	
		double ax = maxBrackets_[i].left;
		double bx = maxBrackets_[i].middle;
		double cx = maxBrackets_[i].right;
		double f1,f2,x0,x1,x2,x3;

		x0 = ax;
		x3 = cx;

		if(fabsp(cx,bx,period_) > fabsp(bx,ax,period_)) {
			x1 = bx;
			x2 = normalize(bx + C*fabsp(cx,bx,period_));
		} else { 
			x2 = bx;
			x1 = normalize(bx - C*fabsp(bx,ax,period_));
		}

		f1 = periodicGaussianMixture(params_,x1,period_);
		f2 = periodicGaussianMixture(params_,x2,period_);

		while(fabsp(x3,x0,period_) > tolerance*(fabs(x1))) {
			if(f2 < f1) {
				double new_x1 = normalize(x2 - C*fabsp(x2,x0,period_));
				shift3(x3,x2,x1,new_x1);
				shift2(f2,f1,periodicGaussianMixture(params_,x1,period_));
			} else {
				double new_x2 = normalize(x1 - C*fabsp(x1,x3,period_));
				shift3(x0,x1,x2,normalize(R*x1+C*x3));
				shift2(f1,f2,periodicGaussianMixture(params_, x2, period_));
			}
		}
		maxima.push_back(x1);
	}
	return maxima;
}

vector<double> MethodsPeriodicGaussian::findMinima() const {
    vector<double> minima;
    vector<double> maxima = findMaxima();
    assert(maxima.size() > 0);
    sort(maxima.begin(), maxima.end());
    const double perturb = 1e-2;
    for(int k=0; k<maxima.size(); k++) {
        double ak = maxima[k];
        double bk = maxima[(k+1)%maxima.size()];
        ak = normalize(ak+perturb);
        bk = normalize(bk-perturb);
        double mk = 0;
        double my = 0;
        int iteration = 0;
        do {
            if(iteration > 1e5) {
                throw(std::runtime_error("Error: findMinima maximized number of iterations reached."));
            }
            else
                iteration++;
            
            //suppose we have the maximas:
            // a b c d, where a<b<c<d 
            // there are minimas in (a,b) (b,c) (c,d) (d,a)
            // the if statement takes care of the (d,a) case
            if(bk < ak)
                mk = normalize((ak+period_+bk)/2.0);
            else
                mk = normalize((ak+bk)/2.0);
            my = periodicGaussianMixtureDx(params_, mk, period_);
            if(my < 0) 
                ak = mk;
            else
                bk = mk;
        } while (fabsp(ak,bk,period_) > 1e-8);
        minima.push_back(mk);
    }
    return minima;
}

}
