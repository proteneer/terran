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

/*
Cluster::Cluster(const vector<vector<double> > &data, const vector<double> &period, const vector<vector<Param> > &initialParams) : 
    dataset_(data),
    period_(period),
    paramset_(initialParams),
    partitions_(paramset_.size()) {

    for(int i=0; i < period_.size(); i++) {
        if(period_[i] < 0) {
            throw(std::runtime_error("period cannot be less than 0"));
        }
    }

    if(data[0].size() != period.size()) {
        throw(std::runtime_error("number of dimensions in data does not match number of dimensions in period!"));
    }
}
*/

Cluster::Cluster(const vector<vector<double> > &data, const vector<double> &period) : 
    dataset_(data),
    period_(period),
    partitions_(period.size()),
    partitioners_(period.size()) {
    if(data.size() == 0) 
        throw(std::runtime_error("Cluster()::Cluster() - input data size cannot 0"));

    cout << data[0].size() << " " << period.size() << endl;


    if(data[0].size() != period.size())
        throw(std::runtime_error("Cluster()::Cluster() - period size does not match data dimension"));


    for(int d=0; d < getNumDimensions(); d++) {
        partitioners_[d] = new PartitionerEM(getDimension(d), period_[d]);

    }
    // set partitioners by default to EM    
}

Cluster::~Cluster() {
    for(int i=0; i < partitioners_.size(); i++) {
        delete partitioners_[i];
        // no need to set to NULL as whole thing will be destroyed.
    }
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

double Cluster::getPeriod(int d) const {
    if(!isPeriodic(d)) {
        stringstream msg;
        msg << "getPeriod() exception, dimension " << d << " is not periodic!";
        throw(std::runtime_error(msg.str()));
    }
    return period_[d];
}

vector<double> Cluster::getPoint(int n) const {
    if(n >= getNumPoints()) {
        throw(std::runtime_error("Cluster::getPoint() out of bounds!"));   
    }
    return dataset_[n];
}

vector<double> Cluster::getDimension(int d) const {
    if(d >= getNumDimensions()) {
        throw(std::runtime_error("Cluster::getMarginalValues() out of bounds!"));   
    }
    vector<double> data(getNumPoints());
    for(int i=0; i<getNumPoints(); i++) {
        data[i] = dataset_[i][d];
    }
    return data;
}

void Cluster::partition(int d) {
    partitions_[d] = partitioners_[d]->partition();
};
/*
vector<Param> Cluster::getParameters(int d) const {
    if(d > getNumDimensions() || d < 0) {
        throw(std::runtime_error("getParameters(), invalid dimension"));
    }
    return paramset_[d];
}
*/

vector<double> Cluster::getPartition(int d) const {
    if(d > getNumDimensions() || d < 0) {
        throw(std::runtime_error("Cluster::getPartitions() - invalid dimension"));
    }
    return partitions_[d];
}

/*
void Cluster::partition(int d, double threshold) {
    vector<Param> params = paramset_[d];
    if(params.size() == 0) {
        throw(std::runtime_error("partition() error: Parameters for dimension 0 have not been optimized!"));
    }
    vector<double> partition;
    if(isPeriodic(d)) {
        const double period = period_[d];
        MethodsPeriodicGaussian mpg(params, period);
        vector<double> minima = mpg.findMinima();
        for(int i=0; i < minima.size(); i++) {
            double val = periodicGaussianMixture(params, minima[i], period);
            if(val < threshold) {
                partition.push_back(minima[i]);
            }
        }
    } else {
        // TODO: implement code to partition non periodic mixture models

    }
    // this is important
    sort(partition.begin(), partition.end());
    partitions_[d] = partition;
}
*/

/*
void Cluster::setParameters(int d, const vector<Param> &p) {
    if(d > getNumDimensions() - 1) {
        throw(std::runtime_error("Dimension out of bounds\n"));
    }
    paramset_[d] = p;
}
*/

void Cluster::setPartition(int d, const vector<double> &p) {
    if(d > getNumDimensions() - 1) {
        throw(std::runtime_error("Dimension out of bounds\n"));
    }
    partitions_[d] = p;
}

vector<int> Cluster::cluster() {
    
    for(int d = 0; d < getNumDimensions(); d++) {
        if(partitions_[d].size() == 0) {
            stringstream msg;
            msg << "Error in cluster(), partitions in dimension " << d << " not set!" << endl;
            throw(std::runtime_error(msg.str()));
        }
    }
    // assign each point to a bucket
    // cannot parallelize easily due to push_back
    map<vector<short>, vector<int> > clusters;
    for(int n = 0; n < getNumPoints(); n++) {
        vector<short> bucket = assign(n);
        cout << n << " ";
        for(int i=0; i < bucket.size(); i++) {
            cout << bucket[i];
        }
        cout << endl;
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

/*
void Cluster::optimizeParameters(int d) {

    // get the marginal distribution for dimension d
    vector<double> data;
    for(int n = 0; n < getNumPoints(); n++) {
        data.push_back(dataset_[n][d]);
    }
    // load initial set of parameters
    const vector<Param> &initialParams = paramset_[d];

    if(isPeriodic(d)) {
        EMPeriodicGaussian epg(data, period_[d]);
        if(initialParams.size() == 0) {
            epg.multiAdaptiveRun(100, 0.1, 0.08, 5, 5);
        } else {
            epg.setParameters(initialParams);
            epg.run();
        }
        paramset_[d] = epg.getParams();
    } else {
        EMGaussian eg(data, initialParams);
        if(initialParams.size() == 0) {
            eg.multiAdaptiveRun(100, 0.1, 0.08, 5, 5);
        } else { 
            eg.setParameters(initialParams);
            eg.run();
        }
        paramset_[d] = eg.getParams();
    }
}
*/

vector<short> Cluster::assign(int pointIndex) const {
    vector<short> bucket(getNumDimensions());
    vector<double> point = dataset_[pointIndex];

    // for each dimension in the point
    for(int d=0; d<getNumDimensions(); d++) {
        // loop over the intervals
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
