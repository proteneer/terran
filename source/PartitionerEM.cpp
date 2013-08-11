#include "PartitionerEM.h"
#include "EMGaussian.h"
#include "EMPeriodicGaussian.h"
#include "MethodsPeriodicGaussian.h"

#include <algorithm>
#include <vector>

using namespace std;
using namespace Terran;

PartitionerEM::PartitionerEM(const vector<double> &dataset, double period) : 
    period_(period),
    em_(NULL),
    MARThreshold_(0.05), 
    MARNumParams_(12),
    MARNumTries_(5),
    partitionCutoff_(0.05) {

    if(period_ != 0) {
        em_ = new EMPeriodicGaussian(dataset, period_);
    } else {
        em_ = new EMGaussian(dataset);
    }

}

PartitionerEM::~PartitionerEM() {
    delete em_;
}

void PartitionerEM::optimizeParameters() {
    em_->multiAdaptiveRun(MARThreshold_, MARNumParams_, MARNumTries_);
}

std::vector<double> PartitionerEM::partition() {
    cout << "optimizing parameters" << endl;
    optimizeParameters();
    cout << "finding low minima points" << endl;
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
            cout << val << " " << partitionCutoff_ << endl;
            if(val < partitionCutoff_) {
                partition.push_back(minima[i]);
            }
        }
    } else {
        // TODO: implement code to partition non periodic mixture models

    }

    // if no minima at all, then just add in a partition randomly for periodic domains

    // for nonperiodic domain, do not add in a partition

    // 
    sort(partition.begin(), partition.end());
    return partition;
}

// this can be used to manipulate the underlying EM object if desired
EM& PartitionerEM::getEM() {
    return *em_;   
};