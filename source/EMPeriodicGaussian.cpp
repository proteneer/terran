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

const int numImages = 10;

void EMPeriodicGaussian::initializePink() {
	// yuck.
	vector<vector<vector<double> > > init(data_.size(), vector<vector<double> >(params_.size(), vector<double>(2*numImages+1, 0)));
	pinkr_ = init;
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
            for(int j=i+1; j != i; j = (j+1) % params_.size()) {
                candidates.push_back(params_[j]);
                Param estimate = estimator(candidates);
                double squaredError = squaredIntegratedError(candidates, estimate);
                if(sqrt(squaredError) < 1e-2) {
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
				pinkr_[n][k][r+numImages] = top/bot;
            }
		}
    }

	for(int n=0; n < data_.size(); n++) {
		double sum = 0;
		for(int k=0; k < params_.size(); k++) {
			for(int r= - numImages; r <= numImages; r++) {
				sum += pinkr_[n][k][r+numImages];   
			}
		}
		if(sum < 0.999) {
			throw std::runtime_error("EMPeriodicGaussian::EStep() failed - pinkr_ no longer sums to 1");
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

/*
void EMPeriodicGaussian::MStep() {
#ifndef NDEBUG
    testIntegrity();
#endif

	cout << "mean" << endl;

    // Compute new mean
    for(int k=0; k<params_.size(); k++) {
        // find two points ak and bk such that
        // ak < bk, and dldu(ak) > dldu(bk)
        // we just divide the domain into 8 equal parts
        // and guess that its bound to cross at one point
        int increment = 16;

        // look for two points that bracket the point
        vector<double> signs(increment);
        vector<double> xvals(increment);
        for(int i=0; i<increment; i++) {
            double x = -period_/2+i*period_/increment;
            signs[i] = dldu(x,k);
            xvals[i] = x;
        }

        bool found = false;
        double ak,bk;
        for(int i=0; i<signs.size(); i++) {
            if(signs[i]>0 && signs[(i+1)%increment]<0) {
                found = true;
                ak = xvals[i];
                bk = xvals[(i+1)%increment];
            }
        }

        if(!found) {
            throw(std::runtime_error("Fatal: Bracket not found!"));
        }

        double mk;
        double my;

        int iteration = 0;
        // bisection with periodic boundaries
        do {
            if(iteration > 1e3) {
                throw(std::runtime_error("EMPeriodicGaussian::MStep() (mean) maximum number of iterations reached."));
            }
             else
                 iteration++;
            if(bk < ak)
                mk = normalize((ak+period_+bk)/2.0);
            else
                mk = normalize((ak+bk)/2.0);
            my = dldu(mk, k);
            if(my > 0)
                ak = mk;
            else
                bk = mk;
        } while(fabs(my) > 1e-5);
        params_[k].u = mk;
    }
	
	cout << "std dev" << endl;

    // Compute new standard deviation and probability
    for(int k=0; k<params_.size(); k++) {
        // bisection root-finding algorithm
        double ak = 0.01;
        double bk = period_;
        double mk;
        double my;

		int count = 0;

		ofstream log("debug.txt");
        do {
			count++;
			if(count > 400) {
				for(double x = 0; x < 10; x++) {
					log << dlds(x, k) << endl;
				}
				throw(std::runtime_error("FATAL"));
			}
            assert(dlds(ak, k) > 0 || isnan(dlds(ak, k)));
            assert(dlds(bk, k) < 0);
            mk = (ak+bk)/2.0;
            my = dlds(mk, k);
            if(my > 0)
                ak = mk;
            else
                bk = mk;
        } while(fabs(my) > 1e-5);
        params_[k].s = mk;

        // Compute new probability
        double denominatorSum = 0;
        for(int n=0; n<data_.size(); n++) {
            denominatorSum += pikn_[n][k];
        }
        params_[k].p = denominatorSum / data_.size();
    }
}
*/

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

/*
double EMPeriodicGaussian::dldu(double uk, int k) const {
    double sk = params_[k].s;
    double s_sum = 0;
    for(int n=0; n<data_.size(); n++) {
        double xn = data_[n];
        double summand = 0;
        double r_sum_rg = 0;
        double r_sum_g = 0;

        double left = uk-7*sk;
        double right = uk+7*sk;
        int lower = floor(left/period_+0.5);
        int upper = floor(right/period_+0.5);

        lower = min(-1, lower);
        upper = max( 1, upper);

        for(int r=lower; r <= upper; r++) {
            r_sum_g  += gaussian(uk+r*period_, sk, xn); 
            r_sum_rg += r*gaussian(uk+r*period_, sk, xn);
        }
        summand = r_sum_rg/r_sum_g;
        if(isnan(summand))
            summand = 0;
        s_sum += pikn_[n][k]*(xn - uk - period_*summand);
    }
    return s_sum;
}

double EMPeriodicGaussian::dlds(double sk, int k) const {
    double s_sum = 0;
    const double uk = params_[k].u;
    for(int n=0; n<data_.size(); n++) {
        const double xn = data_[n];
        double p1 = 0;
        double p2 = 0;
        double p3 = 0;
        p1 = -sk*sk + (xn-uk)*(xn-uk);

        double r_sum_g = 0;
        double r_sum_rg = 0;
        double r_sum_rrg = 0;

        const int D = period_;

        double left = uk-7*sk;
        double right = uk+7*sk;
        int lower = floor(left/period_+0.5);
        int upper = floor(right/period_+0.5);
        lower = min(-1, lower);
        upper = max( 1, upper);
        for(int r=lower; r<=upper; r++) {
            double gauss = gaussian(uk+r*D, sk, xn);
            r_sum_g   += gauss;
            r_sum_rg  += r*gauss;
            r_sum_rrg += r*r*gauss;
        }
        p2 = -2*(xn-uk)*D*r_sum_rg / r_sum_g;
        p3 = D*D*r_sum_rrg / r_sum_g;
        if(isnan(p2))
            p2 = 0;
        if(isnan(p3))
            p3 = 0;
        s_sum += pikn_[n][k]*(p1+p2+p3);
    }
    return s_sum;
}
*/

}
