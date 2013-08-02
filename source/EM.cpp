#include "EM.h"
#include "MathFunctions.h"
#include <iostream>
#include <fstream>


namespace Terran {

using namespace std;

EM::~EM() {

}

EM::EM(const std::vector<double> &data, const std::vector<Param> &params) : 
    data_(data),
    params_(params),
    pikn_(data.size(), std::vector<double>(params.size(),0)) {

    if(data_.size() == 0)
        throw(std::runtime_error("Cannot initialize EM with empty dataset"));

    double psum = 0;

    for(int i=0; i<params_.size(); i++) {
        psum += params_[i].p;
        if(params_[i].p <= 0)
            throw(std::runtime_error("Cannot have p <= 0 in parameters"));
        if(params_[i].p > 1)
            throw(std::runtime_error("Cannot have p > 1 in parameters"));
        if(params_[i].s <= 0)
            throw(std::runtime_error("Cannot have s <= 0 in parameters"));
    }
    if(psum > 1.0001) {
        throw(std::runtime_error("Initial probabilities sum to greater than 1"));
    }
}

void EM::setParams(const std::vector<Param> &input) {
    params_ = input;
}

std::vector<Param> EM::getParams() const {
    return params_;
}

int EM::getDataSize() const {
    return data_.size();
}

double EM::getLikelihood() const {
    double lambda = 0;
    for(int n=0; n<data_.size(); n++) {
        double sum = 0;
        for(int k=0; k<params_.size(); k++) {
            sum += qkn(k,n);
        }
        lambda += log(sum);
    }
    return lambda;
}


static void ppg(const vector<Param> &params, double period, int images) {
    std::ofstream mixture("mixture.dat");
    for(double xn = -period/2; xn < period/2; xn += 0.01) {
        mixture << xn << " " << periodicGaussianMixture(params, xn, period, images) << std::endl;
    }
    std::ofstream mixtureDx("mixtureDx.dat");
    for(double xn = -period/2; xn < period/2; xn += 0.01) {
        mixtureDx << xn << " " << periodicGaussianMixtureDx(params, xn, period, images) << std::endl;
    }
}

bool EM::run(int maxSteps, double tolerance) {
    int steps = 0;
    double likelihood = getLikelihood();
    double likelihoodOld;
    // keep an old copy of params
    vector<Param> paramsOld;

    do {
        likelihoodOld = likelihood;
        paramsOld = params_;
        EStep();
        MStep();   
        steps++;
        likelihood = getLikelihood(); 
        // if the likelihood increased, then we revert back to the old params right before
        // we took the step and break;
        // (the likelihood may increase due to convergence/numerical issues, and is
        // an indication of convergence)

        // likelihood decreases normally if a convergence criterion has been reached
        if(likelihood < likelihoodOld) {
            params_ = paramsOld;
            break; 
        }
    // Stop EM if:
    // a. likelihood reaches the specified tolerance
    // b. maxmimum number of steps reached
    } while(likelihood - likelihoodOld > tolerance && steps < maxSteps);

    return (steps < maxSteps);
}

bool EM::adaptiveRun(int maxSteps, double tolerance, double cutoff) {
    int steps = 0;
    double likelihood = getLikelihood();
    double likelihoodOld;
    // keep an old copy of params
    vector<Param> paramsOld;

    do {
        likelihoodOld = getLikelihood();
        paramsOld = params_;
        EStep();
        MStep();   
        steps++;
        likelihood = getLikelihood(); 
        vector<Param> newParams;
        for(int i=0; i < params_.size(); i++)
            if(params_[i].p > cutoff)
                newParams.push_back(params_[i]);

        if(newParams.size() != params_.size()) {
            params_ = newParams;
        } else if(likelihood < likelihoodOld) {
            params_ = paramsOld;
            break; 
        }
    // Stop EM if:
    // a. likelihood reaches the specified tolerance
    // b. maxmimum number of steps reached
    } while(likelihood-likelihoodOld > tolerance && steps < maxSteps);

    return (steps < maxSteps);
}

/*

move to clustering class

bool EM::multiAdaptiveRun(int maxSteps, double tolerance, double cutoff, int numParams, int numSeeds) {
    vector<Param> bestParams;
    double bestLikelihood;
    int iteration = 0;
    do {
        vector<Param> params;
        for(int i=0; i < numParams; i++) {
            // generate a random number between 0 and 1
            double frac = (double) rand() / (double) RAND_MAX;
            Param p;
            p.p = 1.0/numParams;
            p.u = -period/2 + frac*period; 
            p.s = 0.3;
        }
        params.push_back(p);
          
        iteration++;
    } while(iteration < numSeeds);
}
*/

void EM::EStep() {
    for(int n=0; n<data_.size(); n++) {
        double sum = 0;
        for(int k=0; k<params_.size(); k++) {
            sum += qkn(k,n);
        }
        for(int k=0; k<params_.size(); k++) {
            pikn_[n][k] = qkn(k,n)/sum;
        }
    }
    testIntegrity();
}

void EM::testIntegrity() const {
    double sum = 0;
    for(int k=0; k<params_.size(); k++) {
        for(int n=0; n<data_.size(); n++) {
            sum += pikn_[n][k];   
        }
        if(fabs(sum-1.0) < 1e-7)
            throw(std::runtime_error("pikn no longer sums to 1"));
    }
};

}
