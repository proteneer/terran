#include "PartitionerEM.h"
#include "EMGaussian.h"
#include "EMPeriodicGaussian.h"
#include "MethodsPeriodicGaussian.h"
#include "MethodsGaussian.h"

#include <sstream>
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
	delete em_;
}

void PartitionerEM::optimizeParameters() {
    em_->simpleRun(initialK_);
}

void PartitionerEM::setDataAndPeriod(const vector<double> &data, bool isPeriodic) {
	
	isPeriodic_ = isPeriodic;
	
	if(isPeriodic) {
		for(int i=0; i < data.size(); i++) {
			if( data[i] < -PI || data[i] > PI ) {
				stringstream msg;
				msg << "PartitionerEM::setDataAndPeriod - point " << i << " not in the periodic interval [-PI, PI]";
				throw(std::runtime_error(msg.str()));
			}
		}
	}

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

Partitioner* PartitionerEM::clone(const std::vector<double> &data, bool isPeriodic) {
	PartitionerEM *pem = new PartitionerEM;
	pem->setDataAndPeriod(data, isPeriodic);
	pem->isPeriodic_ = isPeriodic;
	pem->initialK_ = this->initialK_;
	pem->partitionCutoff_ = this->partitionCutoff_;
	return pem;
}

std::vector<double> PartitionerEM::partition() {
	if(em_ == NULL) {
		throw(std::runtime_error("PartitionEM::findLowMinima() - dataset_ has not been initialized"));
	}
    optimizeParameters();
	vector<double> points = findLowMinima();
	return points;
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

static bool compMean(const Param &a, const Param&b) {
	return a.u < b.u;
}

static bool compSig(const Param &a, const Param&b) {
	return a.s < b.s;
}

// get a curve of the model approximating the data
// if periodic, left and right are ignored.
void PartitionerEM::evaluateModel(vector<double> &xvals, vector<double> &yvals, int nsamples) {

	vector<Param> params = em_->getParams();

	double left;
	double right;

	if(isPeriodic_) {
		left = -PI;
		right = PI;
	} else {
		vector<Param>::const_iterator min_mean = min_element(params.begin(), params.end(), compMean);
		vector<Param>::const_iterator max_mean = max_element(params.begin(), params.end(), compMean);
		vector<Param>::const_iterator max_sig  = max_element(params.begin(), params.end(), compSig);

		double left = min_mean->u;
		while(gaussianMixture(params, left) > 1e-2) {
			left -= 3 * max_sig->s;
		}
		double right = max_mean->u;
		while(gaussianMixture(params, right) > 1e-2) {
			right += 3* max_sig->s;
		}

	}
	xvals.resize(0);
	yvals.resize(0);
	for(double x=left; x<right; x += (right-left)/ nsamples) {
		xvals.push_back(x);
		if( isPeriodic_ ) {
			yvals.push_back(periodicGaussianMixture(params, x, 2*PI));
		} else {
			yvals.push_back(gaussianMixture(params, x));
		}
	}
}

// this can be used to manipulate the underlying EM object if desired
EM& PartitionerEM::getEM() {
	if(em_ == NULL) {
		throw(std::runtime_error("PartitionerEM::getEM() - em_ is NULL"));
	}
    return *em_;   
};