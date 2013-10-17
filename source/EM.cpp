#include "EM.h"
#include "MathFunctions.h"
#include <limits>
#include <iostream>
#include <fstream>
#include <sstream>
#include <algorithm>
#include <time.h>

#ifdef _WIN32
#define isinf(x) !_finite(x)
#endif

namespace Terran {

using namespace std;

EM::EM(const std::vector<double> &data) : 
    data_(data),
    //pikn_(data.size(), std::vector<double>(0)),
    maxSteps_(200),
    tolerance_(0.1) {
    if(data_.size() == 0)
        throw(std::runtime_error("Cannot initialize EM with empty dataset"));
}

EM::EM(const std::vector<double> &data, const std::vector<Param> &params) : 
    data_(data),
    params_(params),
    //pikn_(data.size(), std::vector<double>(params.size(),0)),
    maxSteps_(200),
    tolerance_(0.1) {

    if(data_.size() == 0)
        throw(std::runtime_error("Cannot initialize EM with empty dataset"));
    setParameters(params);

}

EM::~EM() {

}

// TODO: make this virtual and initialize pikn
void EM::setParameters(const std::vector<Param> &input) {
    double psum = 0;
    for(int i=0; i<input.size(); i++) {
        psum += input[i].p;
        if(input[i].p <= 0)
            throw(std::runtime_error("Cannot have p <= 0 in parameters"));
        if(input[i].p > 1)
            throw(std::runtime_error("Cannot have p > 1 in parameters"));
        if(input[i].s <= 0)
            throw(std::runtime_error("Cannot have s <= 0 in parameters"));
    }
    if(psum > 1.0000001) {
        throw(std::runtime_error("Initial probabilities sum to greater than 1"));
    }
    params_ = input;
}

std::vector<Param> EM::getParams() const {
    return params_;
}

int EM::getDataSize() const {
    return data_.size();
}

void EM::setMaxSteps(int maxSteps) {
    maxSteps_ = maxSteps;
}

int EM::getMaxSteps() const {
    return maxSteps_;
}

void EM::setTolerance(double tol) {
    tolerance_ = tol;
}

double EM::getTolerance() const {
    return tolerance_;
}

bool EM::run() {
    if(params_.size() == 0) {
        throw(std::runtime_error("EM::run(), parameters are not set"));
    }

	initializePink();

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
    } while(likelihood - likelihoodOld > tolerance_ && steps < maxSteps_);

	destroyPink();

    return (steps < maxSteps_);
}

bool paramSorter(const Param &p1, const Param &p2) {
    return p1.u < p2.u;
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

#include "omp.h"

bool EM::simpleRun(unsigned int numParams) {
    
	if(numParams > data_.size()) {
        throw(std::runtime_error("EM::simpleRun(), numParams > number of data points"));
    }

    // initialize parameters by sampling from the data
    vector<double> randomSample = data_;
    random_shuffle(randomSample.begin(), randomSample.end());
    randomSample.resize(numParams);
    params_.resize(numParams);

    for(int i=0; i < randomSample.size(); i++) {
        params_[i].p = (double) 1 / randomSample.size();
        params_[i].u = randomSample[i];
        params_[i].s = 0.1*domainLength();
    }

	initializePink();
	int steps = 0;
    double likelihood = getLikelihood();
    double likelihoodOld;
    do {
        vector<Param> paramsOld = params_;
		likelihoodOld = likelihood;
		EStep();
		MStep();
        steps++;
        likelihood = getLikelihood(); 
        if(steps >= maxSteps_) {
            break;
        }

        int initialSize = params_.size();
		mergeParams();
        if(initialSize != params_.size()) {
            continue;
        }
        // rethink termination criteria if there are merges happening
        if(fabs(likelihood - likelihoodOld) < tolerance_) {
            if(likelihood < likelihoodOld) {
                params_ = paramsOld;
            }
            break;
        }

    // Stop EM if:
    // a. likelihood reaches the specified tolerance
    // b. maximimum number of steps reached
    } while(true);
	
	destroyPink();

    return steps < maxSteps_;
}

}
