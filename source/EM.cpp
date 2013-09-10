#include "EM.h"
#include "MathFunctions.h"
#include <limits>
#include <iostream>
#include <fstream>
#include <sstream>
#include <algorithm>


namespace Terran {

using namespace std;

EM::~EM() {

}

EM::EM(const std::vector<double> &data) : 
    data_(data),
    pikn_(data.size(), std::vector<double>(0)),
    maxSteps_(100),
    tolerance_(0.1) {
    if(data_.size() == 0)
        throw(std::runtime_error("Cannot initialize EM with empty dataset"));
}

EM::EM(const std::vector<double> &data, const std::vector<Param> &params) : 
    data_(data),
    params_(params),
    pikn_(data.size(), std::vector<double>(params.size(),0)),
    maxSteps_(100),
    tolerance_(0.1) {

    if(data_.size() == 0)
        throw(std::runtime_error("Cannot initialize EM with empty dataset"));

    setParameters(params);
}

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
    vector<vector<double> > temp(data_.size(), std::vector<double>(params_.size(),0));
    pikn_ = temp; 
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

    return (steps < maxSteps_);
}

// the number of components decrease. 
bool EM::adaptiveRun(double cutoff) {

    if(params_.size() == 0) {
        throw(std::runtime_error("EM::adaptiveRun(), parameters are not set"));
    }

	for(int i=0; i < data_.size(); i++) {
		pikn_[i].resize(params_.size(), 0);
	}

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
    // b. maximimum number of steps reached
    } while(likelihood-likelihoodOld > tolerance_ && steps < maxSteps_);

    return (steps < maxSteps_);
}

// TODO: set pikn size properly!! 
void EM::multiAdaptiveRun(double cutoff, int numParams, int numTries) {
    vector<Param> bestParams;

    double bestLikelihood = -numeric_limits<double>::max();
    int attempts = 0;

    do {
        // initialize a set of random parameters
        vector<Param> params;
        vector<double> mean = sampleDomain(numParams);
        for(int i=0; i < mean.size(); i++) {
            Param p;
            p.p = 1.0/numParams;
            p.u = mean[i];
            // todo, should s be variable as well?
            p.s = 0.3;
            params.push_back(p);
        }
        setParameters(params);
        adaptiveRun(cutoff);
        double newLikelihood = getLikelihood();

        if(newLikelihood > bestLikelihood) {
            bestLikelihood = newLikelihood;
            bestParams = getParams(); 
        }
        /*
        vector<Param> testParams = getParams(); 
        for(int i=0; i < testParams.size(); i++) {
            cout << testParams[i].p << " " << testParams[i].u << " " << testParams[i].s << endl;
        }
        cout << "Likelihood: " << getLikelihood() << endl;
        */
        attempts++;
    } while(attempts < numTries);
    params_ = bestParams;
}

void EM::kernelAdaptiveRun() {

	cout << "starting KAR" << endl;

	// bandwidth estimation needs to be subclassed and computed
	// periodic and nonperiodic spaces compute differently presumably
	double bandwidth = 0.3;

	vector<double> randsample = data_;
    srand(time(0));
	random_shuffle(randsample.begin(), randsample.end());
	unsigned int numParams = min((int)randsample.size(), 1000);

	params_.resize(numParams);
	for(int k=0; k < numParams; k++) {
		Param p;
		p.u = data_[k];
		p.s = bandwidth;
		p.p = (double) 1 / (double) params_.size();
		params_[k] = p;
	}

	for(int i=0; i < data_.size(); i++) {
		pikn_[i].resize(params_.size(), 0);
	}

	// adaptive run using a varying cutoff (slowly increasing);
	// invariant: cutoff must always be > 1 / num_components
	// cutoff set cutoff to 1/10th of 1/num_components?

	int steps = 0;

    double likelihood = getLikelihood();
    double likelihoodOld;

    // keep an old copy of params
    vector<Param> paramsOld;
/*
	    for(int i=0; i < params_.size(); i++) {
			cout << params_[i].p << " " << params_[i].u << " " << params_[i].s <<  endl;
		}
		return;
*/
    do {
		// slowly increase params size as needed

        //cout << "Step: " << steps << " " << params_.size() << endl;

		double cutoff = min(0.5/params_.size(),0.02);
        likelihoodOld = getLikelihood();
        paramsOld = params_;
		
        //cout << "e" << endl;
        EStep();
        //for(int i=0; i < params_.size(); i++) {
        //    cout << params_[i].p << " " << params_[i].u << " " << params_[i].s << endl;
        //}
		//cout << "m" << endl;
        MStep(); 
        //for(int i=0; i < params_.size(); i++) {
        //    cout << params_[i].p << " " << params_[i].u << " " << params_[i].s << endl;
        //}
        
        steps++;
        likelihood = getLikelihood(); 
        vector<Param> newParams;
        double min = 999;
        for(int i=0; i < params_.size(); i++) {
			//cout << params_[i].p << endl;

            if(params_[i].p < min)
                min = params_[i].p;

            if(params_[i].p > cutoff)
                newParams.push_back(params_[i]);
		}

        if(newParams.size() != params_.size()) {
			//cout << "deleted " << params_.size() - newParams.size() << " components." << endl;
            params_ = newParams;
        } else if(likelihood < likelihoodOld) {
            params_ = paramsOld;
            break; 
        }

        //cout << likelihood << " " << likelihoodOld << " " << likelihood-likelihoodOld << endl;
    // Stop EM if:
    // a. likelihood reaches the specified tolerance
    // b. maximimum number of steps reached
    } while( isinf(likelihood) ? true : (likelihood-likelihoodOld > tolerance_ && steps < maxSteps_));


}

void EM::EStep() {
    for(int n=0; n<data_.size(); n++) {
        double sum = 0;
        for(int k=0; k<params_.size(); k++) {
            sum += qkn(k,n);
        }

        for(int k=0; k<params_.size(); k++) {
            if(sum > 1e-7)
                pikn_[n][k] = qkn(k,n)/sum;
            else
                pikn_[n][k] = 0;
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
