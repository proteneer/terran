#include "PartitionerEM.h"
#include "EMGaussian.h"
#include "EMPeriodicGaussian.h"
#include "MethodsPeriodicGaussian.h"
#include "MethodsGaussian.h"

#include <algorithm>

using namespace std;
using namespace Terran;

PartitionerEM::PartitionerEM() :
	Partitioner(),
	partitionCutoff_(0.05),
	em_(NULL),
	initialK_(25) {

}

PartitionerEM::~PartitionerEM() {
    if(em_ != NULL)
		delete em_;
}

void PartitionerEM::optimizeParameters() {
    em_->simpleRun(initialK_);
}

void PartitionerEM::setDataAndPeriod(const vector<double> &data, bool isPeriodic) {
	
	isPeriodic_ = isPeriodic;
	
	// delete the old em_ object
	if(em_ != NULL) {
		delete em_;
		em_ = NULL;
	}

	// instantiate a new em_ object
    if(isPeriodic) {
        em_ = new EMPeriodicGaussian(data, 2*PI);
    } else {
        em_ = new EMGaussian(data);
    }

}

std::vector<double> PartitionerEM::partition() {
	if(em_ == NULL) {
		throw(std::runtime_error("PartitionEM::findLowMinima() - dataset_ has not been initialized"));
	}
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

double PartitionerEM::getPartitionCutoff() const {
	return partitionCutoff_;
}

void PartitionerEM::setInitialK(int count) {
	initialK_ = count;
}

int PartitionerEM::getInitialK() const {
	return initialK_;
}

// this can be used to manipulate the underlying EM object if desired
EM& PartitionerEM::getEM() {
	if(em_ == NULL) {
		throw(std::runtime_error("PartitionerEM::getEM() - em_ is NULL"));
	}
    return *em_;   
};