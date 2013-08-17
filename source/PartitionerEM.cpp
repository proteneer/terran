#include "PartitionerEM.h"
#include "EMGaussian.h"
#include "EMPeriodicGaussian.h"
#include "MethodsPeriodicGaussian.h"
#include "MethodsGaussian.h"

#include <algorithm>

using namespace std;
using namespace Terran;

PartitionerEM::PartitionerEM(const vector<double> &dataset, double period) : 
    Partitioner(dataset, period),
    MARThreshold_(0.05), 
    MARNumParams_(12),
    MARNumTries_(5),
    partitionCutoff_(0.05),
    em_(NULL) {

    if(isPeriodic()) {
        em_ = new EMPeriodicGaussian(dataset_, period_);
    } else {
        em_ = new EMGaussian(dataset_);
    }

}

PartitionerEM::~PartitionerEM() {
    delete em_;
}

void PartitionerEM::setMARThreshold(double val) {
    MARThreshold_ = val;
}

void PartitionerEM::setMARNumParams(int val) {
    MARNumParams_ = val;
}

void PartitionerEM::setMARNumTries(int val) {
    MARNumTries_ = val;
}

void PartitionerEM::setPartitionCutoff(double val) {
    partitionCutoff_ = val;
}

double PartitionerEM::getMARThreshold() const {
    return MARThreshold_;
}

int PartitionerEM::getMARNumParams() const {
    return MARNumParams_;
}

int PartitionerEM::getMARNumTries() const {
    return MARNumTries_;
}

double PartitionerEM::getPartitionCutoff() const {
    return partitionCutoff_;
}

void PartitionerEM::optimizeParameters() {
    em_->multiAdaptiveRun(MARThreshold_, MARNumParams_, MARNumTries_);
}

std::vector<double> PartitionerEM::partition() {
    optimizeParameters();
    return findLowMinima();
}

bool PartitionerEM::isPeriodic() const {
    if(period_ == 0) {
        return false;
    } else {
        return true;
    }
}

vector<double> PartitionerEM::findLowMinima() const {
    vector<Param> params = em_->getParams();
    if(params.size() == 0) {
        throw(std::runtime_error("PartitionEM::partition() - Parameters have not been optimized!"));
    }

    vector<double> partition;
    if(isPeriodic()) {
        const double period = period_;
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

// this can be used to manipulate the underlying EM object if desired
EM& PartitionerEM::getEM() {
    return *em_;   
};