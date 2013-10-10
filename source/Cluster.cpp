#include <stdexcept>
#include <sstream>
#include <map>
#include <algorithm>

#include "Cluster.h"
#include "MethodsPeriodicGaussian.h"
#include "EMPeriodicGaussian.h"
#include "EMGaussian.h"

using namespace std;

namespace Terran {

Cluster::Cluster(const vector<vector<double> > &data, const vector<double> &period) : 
    dataset_(data),
    period_(period),
    partitions_(period.size()),
    partitionMethod_("EM"),
    subsampleCount_(min(3000,(int)data.size())) {

    if(data.size() == 0) 
        throw(std::runtime_error("Cluster()::Cluster() - input data size cannot be 0"));

    if(data[0].size() != period.size())
        throw(std::runtime_error("Cluster()::Cluster() - period size does not match data dimension"));

}

Cluster::~Cluster() {

}

int Cluster::getNumDimensions() const {
    return dataset_[0].size();
}

int Cluster::getNumPoints() const {
    return dataset_.size();
}

bool Cluster::isPeriodic(int d) const {
    if(period_[d] == 0)
        return false;
    else
        return true;
}

void Cluster::setSubsampleCount(int count) {
    subsampleCount_ = count;
}

int Cluster::getSubsampleCount() const {
    return subsampleCount_;
}

double Cluster::getPeriod(int d) const {
    if(!isPeriodic(d)) {
        stringstream msg;
        msg << "Cluster::getPeriod() - dimension " << d << " is not periodic!";
        throw(std::runtime_error(msg.str()));
    }
    return period_[d];
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

void Cluster::setPartitionMethod(string type) {
    if(type.compare("EM") == 0) {
        partitionMethod_ = type;
    } else {
        throw std::runtime_error("Cluster::setPartitionMethod() - unknown partition method");
    }
}

void Cluster::partition(int d) {
    if(partitionMethod_.compare("EM") == 0) {
        PartitionerEM pem(getDimension(d), period_[d]);
        partitions_[d] = pem.partition();
    }
};

vector<int> Cluster::cluster() {
    // assign each point to a bucket
    map<vector<short>, vector<int> > clusters;
    for(int n = 0; n < getNumPoints(); n++) {
        vector<short> bucket = assign(n);
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

vector<short> Cluster::assign(int pointIndex) const {
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
