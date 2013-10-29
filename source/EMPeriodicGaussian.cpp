#include "EMPeriodicGaussian.h"
#include <math.h>
#include <assert.h>
#include <algorithm>
#include <iostream>
#include <fstream>
#include <sstream>
#include <complex>

using namespace std;

namespace Terran {

EMPeriodicGaussian::EMPeriodicGaussian(const vector<double> &data,  double period) : 
    EM(data), 
    period_(period) {
    for(int i=0; i<params_.size(); i++) {
        if(params_[i].s > period_)
            throw(std::runtime_error("Cannot have s > period in parameters"));
        if(params_[i].u < -period_/2)
            throw(std::runtime_error("Cannot have u < -period/2"));
        if(params_[i].u > period_/2)
            throw(std::runtime_error("Cannot have u > period/2"));
    }
	initializePink();
}

EMPeriodicGaussian::EMPeriodicGaussian(const vector<double> &data, const vector<Param> &params, double period) : 
    EM(data, params), 
    period_(period) {
    for(int i=0; i<params.size(); i++) {
        if(params_[i].s > period_)
            throw(std::runtime_error("Cannot have s > period in parameters"));
        if(params_[i].u < -period_/2)
            throw(std::runtime_error("Cannot have u < -period/2"));
        if(params_[i].u > period_/2)
            throw(std::runtime_error("Cannot have u > period/2"));
    }
	initializePink();
}

EMPeriodicGaussian::~EMPeriodicGaussian() {

}

// estimator for periodic gaussian
static Param estimator(const vector<Param> &source) {
    vector<Param> params(source);
    Param estimate;
    estimate.p = 0;
    estimate.u = 0;
    estimate.s = 0;

    // calculate new p
    for(int i=0; i < params.size(); i++) {
        estimate.p += params[i].p;
    }

    // renormalize the weights based on copy
    for(int i=0; i < params.size(); i++) {
        params[i].p = params[i].p/estimate.p;
    }

	// calculate new u and s using phase information
	double real = 0;
	double imag = 0;
	for(int i=0; i < params.size(); i++) {
        real += params[i].p*exp(-(params[i].s*params[i].s)/2)*cos(params[i].u);
		imag += params[i].p*exp(-(params[i].s*params[i].s)/2)*sin(params[i].u);
    }
	complex<double> z(real, imag);
	estimate.u = arg(z);
	double R = abs(z);
	estimate.s = sqrt(log(1/(R*R)));
	return estimate;
}

static double normalCDF(double x, double u, double s) {
	double prefix = 0.5;
	double suffix = 1.0+erf((x-u)/(s*sqrt(2.0)));
    return prefix*suffix;
}

static double normalizer(const Param &a, const Param &b) {
    double p1 = a.p;
    double p2 = b.p;
    double u1 = a.u;
    double u2 = b.u;
    double s1 = a.s;
    double s2 = b.s;
    double sum = 0;
    for(int r1=-6; r1 <=6; r1++) {
        for(int r2=-6; r2 <=6; r2++) {
            double u1_n = u1+r1*2*PI;
            double u2_n = u2+r2*2*PI;
            double top = exp(-0.5*((u1_n-u2_n)*(u1_n-u2_n))/(s1*s1+s2*s2));
            double bot = sqrt(2*PI*(s1*s1+s2*s2));
            double u_new = (u1_n*s2*s2+u2_n*s1*s1)/(s2*s2+s1*s1);
            double s_new = (s1*s2)/sqrt(s1*s1+s2*s2);
            sum += (top/bot) * (normalCDF(PI,u_new,s_new)-normalCDF(-PI,u_new,s_new));
        }
    }
    return p1*p2*sum;
}

// integrated squared error 
static double squaredIntegratedError(const vector<Param> &params, const Param &estimate) {
    double left = 0;
    double middle = 0;
    double right = 0;
    for(int i=0; i < params.size(); i++) {
        double p1 = params[i].p;
        double u1 = params[i].u;
        double s1 = params[i].s;
        for(int j=0; j < params.size(); j++) {
            double p2 = params[j].p;
            double u2 = params[j].u;
            double s2 = params[j].s;
            left += normalizer(params[i], params[j]);
        }
        middle += normalizer(params[i], estimate);
    }
    right = normalizer(estimate, estimate);

    return left - 2*middle + right;
}

static bool paramComparator(const Param &a, const Param&b) {
    return a.u < b.u;
}

const int numImages = 7;

void EMPeriodicGaussian::initializePink() {
	// yuck.
	vector<vector<vector<double> > > init(data_.size(), vector<vector<double> >(params_.size(), vector<double>(2*numImages+1, 0)));
	pinkr_ = init;
}

void EMPeriodicGaussian::destroyPink() {
	pinkr_.resize(0);
}

void EMPeriodicGaussian::mergeParams() {

    sort(params_.begin(), params_.end(), paramComparator);
    vector<bool> skip(params_.size(), 0);
    // final set of parameters
    vector<Param> refined;
    for(int i=0; i < params_.size(); i++) {
        // if this parameter has not already been merged
        if(!skip[i]) {
            // find longest continuous sequence of parameters 
            // that can be merged into a single parameter
            bool hasMergedOnce = false;
            vector<Param> candidates;
            Param best = params_[i]; 
            candidates.push_back(params_[i]);
            for(int j=(i+1) % params_.size(); j != i; j = (j+1) % params_.size()) {
                candidates.push_back(params_[j]);
                Param estimate = estimator(candidates);
                double squaredError = squaredIntegratedError(candidates, estimate);
                if(sqrt(squaredError) < 5e-3) {
                    hasMergedOnce = true;
                    best = estimate;
                    // inclusive merge of all continuous indices
                    for(int k=i; k<=j; k++) {
                        skip[k] = true;
                    }
                } else {
                    if(hasMergedOnce) {
                        break;
                    }
                }

            }
            refined.push_back(best);
        }
    }
     
    if(refined.size() != params_.size()) {
        params_ = refined;
    }
}

void EMPeriodicGaussian::EStep() {
	for(int n=0 ; n < data_.size(); n++) {
		double bot = periodicGaussianMixture(params_, data_[n], period_);
		for(int k=0; k < params_.size(); k++) {
			double topsum = 0;
            for(int r = -numImages; r <= numImages; r++) {
				double top = params_[k].p*gaussian(params_[k].u, params_[k].s, data_[n]-period_*r);
				//pinkr_[n][k][r+numImages] = top/bot;
                if(bot > 1e-7)
					pinkr_[n][k][r+numImages] = top/bot;
				else
					pinkr_[n][k][r+numImages] = 0;
			
			}
		}
    }	
}

void EMPeriodicGaussian::MStep() {
    for(int k=0; k < params_.size(); k++) {
		// compute new probability and normalization constant
		double normalization;
        double sum = 0;
        for(int n=0; n < data_.size(); n++) {   
            for(int r = -numImages; r <= numImages; r++) {
                sum += pinkr_[n][k][r+numImages];
            }
        }

        normalization = sum;
		params_[k].p = sum/data_.size();

		// compute new mean
        double numerator = 0;
        for(int n=0; n < data_.size(); n++) {
            for(int r=-numImages; r <= numImages; r++) {
                numerator += pinkr_[n][k][r+numImages]*(data_[n]-r*period_);
            }
        }
		params_[k].u = numerator/normalization;

		// compute new standard deviation
        numerator = 0;
        for(int n=0; n < data_.size(); n++) {
            for(int r=-numImages; r <= numImages; r++) {
                double a = data_[n]-params_[k].u-r*period_;
                numerator += pinkr_[n][k][r+numImages]*(a*a);
            }
        }
        params_[k].s = sqrt(numerator/normalization);
	}
}

double EMPeriodicGaussian::domainLength() const {
    return period_;
}

double EMPeriodicGaussian::qkn(int k, int n) const {
    double xn = data_[n];
    double pk = params_[k].p;
    double uk = params_[k].u;
    double sk = params_[k].s;
    return pk*periodicGaussian(uk,sk,xn,period_);
}

}
