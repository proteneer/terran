#include "EMGaussian.h"
#include "MathFunctions.h"

#include <algorithm>
#include <limits>

namespace Terran {

EMGaussian::EMGaussian(const std::vector<double> &data) : 
    EM(data),
	pink_(data.size(), std::vector<double>(0)) {

}

EMGaussian::EMGaussian(const std::vector<double> &data, const std::vector<Param> &params) : 
    EM(data, params),
	pink_(data.size(), std::vector<double>(0)) {

}

EMGaussian::~EMGaussian() { 

}

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

	// calculate new u
    for(int i=0; i < params.size(); i++) {
        estimate.u += params[i].p*params[i].u;
    }
    // calculate new s
    for(int i=0; i < params.size(); i++) {
        double p = params[i].p;
        double u = params[i].u;
        double s = params[i].s;
        double u_new = estimate.u;
        estimate.s += p*(s*s+u*u-2*u_new*u+u_new*u_new);
    }
    estimate.s = sqrt(estimate.s);
    return estimate;
}

static double normalizer(Param a, Param b) {
    double prefix = a.p * b.p;
    double top = exp(-0.5*((a.u-b.u)*(a.u-b.u))/(a.s*a.s+b.s*b.s));
    double bot = sqrt(2*PI*(a.s*a.s+b.s*b.s));
    return prefix*top/bot;
}

// closed form expression of the integrated squared error of the parameters and the estimate
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

void EMGaussian::mergeParams() {

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
            for(int j=i+1; j < params_.size(); j++) {
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

void EMGaussian::initializePink() {
	for(int i=0; i < pink_.size(); i++) {
        pink_[i].resize(params_.size(), 0);
    }
}

void EMGaussian::destroyPink() {
	pink_.resize(0);
}

void EMGaussian::EStep() {
    for(int n=0; n<data_.size(); n++) {
        double denominator = gaussianMixture(params_, data_[n]);
        for(int k=0; k<params_.size(); k++) {
            if(denominator > 1e-7)
                pink_[n][k] = qkn(k,n)/denominator;
            else
                pink_[n][k] = 0;
        }
    }
}

void EMGaussian::MStep() {
    for(int k=0; k<params_.size(); k++) {
		//cout << k << endl;
        // Compute new mean
        double numeratorSum = 0;
        double denominatorSum = 0;
        for(int n=0; n<data_.size(); n++) {
            numeratorSum += pink_[n][k]*data_[n];
        }
        for(int n=0; n<data_.size(); n++) {
            denominatorSum += pink_[n][k];
        }
        params_[k].u = numeratorSum / denominatorSum;

        // Compute new standard deviation
        numeratorSum = 0;
        for(int n=0; n<data_.size(); n++) {
            double dx = data_[n]-params_[k].u;
            numeratorSum += pink_[n][k]*dx*dx;
        }
        params_[k].s = sqrt(numeratorSum / denominatorSum);

        // Compute new probability
        params_[k].p = denominatorSum / data_.size(); 
    }
}

double EMGaussian::qkn(int k, int n) const {
   double pk = params_[k].p;
   double uk = params_[k].u;
   double sk = params_[k].s;    
   double xn = data_[n];
   return pk * gaussian(uk, sk, xn);
}

double EMGaussian::domainLength() const {
    double min =  numeric_limits<double>::max();
    double max = -numeric_limits<double>::max();
    for(int i=0; i < data_.size(); i++) {
        min = (data_[i] < min) ? data_[i] : min;
        max = (data_[i] > max) ? data_[i] : max;
    }
    return max-min;
}

}
