#include "EM.h"
#include "MathFunctions.h"

EM::~EM() {

}

EM::EM(const std::vector<double> &data, const std::vector<Param> &params) : 
    data_(data),
    params_(params),
    pikn_(data.size(), std::vector<double>(params.size(),0)) {

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
    if(lambda < 0) {
        throw(std::runtime_error("Negative likelihood found! Aborting...\n"));
    }
    return lambda;
}

bool EM::run(int maxSteps, double tolerance) {
    int steps = 0;
    double likelihood = getLikelihood();
    double likelihoodOld = likelihood;
    do {
        likelihoodOld = likelihood;
        EStep();
        MStep();   
        steps++;
        likelihood = getLikelihood(); 
    } while(fabs(likelihoodOld - likelihood) > tolerance && steps < maxSteps);
    return (steps < maxSteps);
}


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
        
        if(fabs(sum-1.0) < 0.01) {
            throw(std::runtime_error("pikn no longer sums to 1"));
        }
    }
};
