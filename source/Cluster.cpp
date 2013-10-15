#include <stdexcept>
#include <sstream>
#include <map>
#include <algorithm>

#include "Cluster.h"
#include "MethodsPeriodicGaussian.h"
#include "EMPeriodicGaussian.h"
#include "EMGaussian.h"
#include "PartitionerEM.h"

using namespace std;

namespace Terran {

Cluster::Cluster(const vector<vector<double> > &data, const vector<int> &period) : 
    dataset_(data),
    period_(period),
    partitions_(period.size()),
	partitionFlag_(period.size(), 0),
    subsampleCount_(min(3000,(int)data.size())) {

	partitioner_ = new PartitionerEM();

	initialize();

}

Cluster::Cluster(const vector<vector<double> > &data, const vector<int> &period, Partitioner* partitioner) :
    dataset_(data),
    period_(period),
    partitions_(period.size()),
	partitionFlag_(period.size(), 0),
    subsampleCount_(min(3000,(int)data.size())),
	partitioner_(partitioner) {
	
	initialize();

}

void Cluster::initialize() {
	
	for(int i=0; i < period_.size(); i++) {
		if(period_[i] != 1 && period_[i] != 0) {
			throw(std::runtime_error("period must either be zero (false), or one (true)"));
		}
	}

	for(int n=1; n < dataset_.size(); n++) {
		if(dataset_[n].size() != dataset_[n-1].size()) {
			throw(std::runtime_error("Cluster::Cluster() - not all points have the same dimensions"));
		}
	}

	for(int n=0; n < dataset_.size(); n++) {
		for(int d=0; d < dataset_[n].size(); d++) {
			if(period_[d]) {
				if(dataset_[n][d] < -PI || dataset_[n][d] > PI) {
					stringstream error;
					error << "Cluster::Cluster() - dimension " << d << " is periodic, but the angles are not in the range [-PI, to PI]" << endl;
					throw(std::runtime_error(error.str()));
				}
			}
		}
	}
		
    if(dataset_.size() == 0) 
        throw(std::runtime_error("Cluster()::Cluster() - input data size cannot be 0"));

    if(dataset_[0].size() != period_.size())
        throw(std::runtime_error("Cluster()::Cluster() - period size does not match data dimension"));
}

Cluster::~Cluster() {
	delete partitioner_;
}

int Cluster::getNumDimensions() const {
    return dataset_[0].size();
}

int Cluster::getNumPoints() const {
    return dataset_.size();
}

bool Cluster::isPeriodic(int d) const {
	return period_[d];
}

void Cluster::setSubsampleCount(int count) {
    if(count < getNumPoints()) {
		throw(std::runtime_error("Cluster::setSubsampleCount() - subsample count can not be less than number of points")); 
	}
	subsampleCount_ = count;
}

int Cluster::getSubsampleCount() const {
    return subsampleCount_;
}

vector<double> Cluster::getPoint(int n) const {
    if(n >= getNumPoints()) {
        throw(std::runtime_error("Cluster::getPoint() - n out of bounds!"));   
    }
    return dataset_[n];
}

vector<double> Cluster::getDimension(int d) const {
    if(d >= getNumDimensions()) {
        throw(std::runtime_error("Cluster::getDimension() - d out of bounds!"));   
    }
    vector<double> data(getNumPoints());
    for(int i=0; i<getNumPoints(); i++) {
        data[i] = dataset_[i][d];
    }
    random_shuffle(data.begin(), data.end());
    data.resize(subsampleCount_);
    return data;
}

vector<double> Cluster::getPartition(int d) const {
    if(d > getNumDimensions() || d < 0) {
        throw(std::runtime_error("Cluster::getPartitions() - invalid dimension"));
    }
    return partitions_[d];
}

void Cluster::setPartition(int d, const vector<double> &p) {
    if(d > getNumDimensions() - 1) {
        throw(std::runtime_error("Dimension out of bounds\n"));
    }
    partitions_[d] = p;
}

void Cluster::partition(int d) {
	partitioner_->setDataAndPeriod(getDimension(d), period_[d]);
	partitions_[d] = partitioner_->partition();
	partitionFlag_[d] = true;
};

void Cluster::partitionAll() {
	for(int d=0; d < getNumDimensions(); d++) {
		partition(d);
	}
}

vector<int> Cluster::assign() {

	for(int d=0; d < getNumDimensions(); d++) {
		if(partitionFlag_[d] == false) {
			stringstream errmsg;
			errmsg << "Cluster::assign() - dimension " << d << " has not been partitioned yet!";
			throw(std::runtime_error(errmsg.str()));
		}
	}

    // assign each point to a bucket
    map<vector<short>, vector<int> > clusters;
    for(int n = 0; n < getNumPoints(); n++) {
        vector<short> bucket = findBucket(n);
        clusters[bucket].push_back(n);
    }

    // loop over the buckets and set to cluster
    int clusterIndex = 0;
    vector<int> assignment(getNumPoints(),-1);
    for(map<vector<short>, vector<int> >::const_iterator it = clusters.begin();
        it != clusters.end(); it++) {
        vector<int> pointIndices = it->second;
        for(int j=0; j<pointIndices.size(); j++) {
            assignment[pointIndices[j]] = clusterIndex;
        }
        clusterIndex++;
    }
    return assignment;
}

vector<short> Cluster::findBucket(int pointIndex) const {
    // the initialization to zero for each dimension in the bucket is important
    // as it takes care of the case when no partitions exist for that dimension
    vector<short> bucket(getNumDimensions(), 0);
    vector<double> point = dataset_[pointIndex];

    // for each dimension in the point
    for(int d=0; d<getNumDimensions(); d++) {
        for(int j=0; j<partitions_[d].size(); j++) {
            if(point[d] < partitions_[d][j]) {
                bucket[d] = j;
                break;
            }
            if(isPeriodic(d))
                bucket[d] = 0;
            else
                bucket[d] = j+1;
        }
    }
    return bucket;
}

} // namespace Terran
