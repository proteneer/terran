#include "EMPeriodicGaussian.h"
#include <math.h>
#include <assert.h>
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
        real += params[i].p*exp(-(params[i].s*params[i].s)*cos(params[i].u));
		imag += params[i].p*exp(-(params[i].s*params[i].s)*sin(params[i].u));
    }
	complex<double> z(real, imag);
	estimate.u = arg(z);
	double R = abs(z);
	estimate.s = sqrt(log(1/(R*R)));

	return estimate;
}

static double normalCDF(double x, double u, double s) {
	double prefix = 0.5;
	double suffix = 1.0+erf((x-u)/(s*sqrt(2)));
}

static double normalizer(double p1, double u1, double s1, double p2, double u2, double s2) {
	double sum = 0
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

void EMPeriodicGaussian::MStep() {
#ifndef NDEBUG
    testIntegrity();
#endif

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

    // Compute new standard deviation and probability
    for(int k=0; k<params_.size(); k++) {
        // bisection root-finding algorithm
        double ak = 0.01;
        double bk = period_;
        double mk;
        double my;
        do {
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
        /*
        int images_ = 50;
        for(int r=-images_; r<=images_; r++) {
            r_sum_g  += gaussian(uk+r*period_, sk, xn); 
            r_sum_rg += r*gaussian(uk+r*period_, sk, xn);
        }*/
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

}
