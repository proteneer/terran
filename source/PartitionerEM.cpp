#include "PartitionerEM.h"
#include "EMGaussian.h"
#include "EMPeriodicGaussian.h"
#include "MethodsPeriodicGaussian.h"
#include "MethodsGaussian.h"

#include <algorithm>

using namespace std;
using namespace Terran;

PartitionerEM::PartitionerEM(const vector<double> &dataset, bool isPeriodic) : 
    Partitioner(dataset, isPeriodic),
    partitionCutoff_(0.05),
    em_(NULL) {

    if(isPeriodic_) {
        em_ = new EMPeriodicGaussian(dataset_, 2*PI);
    } else {
        em_ = new EMGaussian(dataset_);
    }

}

PartitionerEM::~PartitionerEM() {
    delete em_;
}

void PartitionerEM::optimizeParameters() {
    em_->simpleRun(25);
}

std::vector<double> PartitionerEM::partition() {
    optimizeParameters();
    return findLowMinima();
}

vector<double> PartitionerEM::findLowMinima() const {
    vector<Param> params = em_->getParams();
    if(params.size() == 0) {
        throw(std::runtime_error("PartitionEM::findLowMinima() - Parameters do not exist!"));
    }

    vector<double> partition;
    if(isPeriodic_) {
        const double period = 2*PI;
        MethodsPeriodicGaussian mpg(params, period);
        vector<double> minima = mpg.findMinima();
        for(int i=0; i < minima.size(); i++) {
            double val = periodicGaussianMixture(params, minima[i], period);
            if(val < partitionCutoff_) {
                partition.push_back(minima[i]);
            }
        }
    } else {
        MethodsGaussian mg(params);
        vector<double> minima = mg.findMinima();
        for(int i=0; i < minima.size(); i++) {
            double val = gaussianMixture(params, minima[i]);
            if(val < partitionCutoff_) {
                partition.push_back(minima[i]);
            }
        }
    }

    sort(partition.begin(), partition.end());
    return partition;
}

void PartitionerEM::setPartitionCutoff(double cutoff) {
    partitionCutoff_ = cutoff;
}

// this can be used to manipulate the underlying EM object if desired
EM& PartitionerEM::getEM() {
    return *em_;   
};