#include <stdexcept>
#include <sstream>
#include <map>

#include "Cluster.h"
#include "MethodsPeriodicGaussian.h"
#include "EMPeriodicGaussian.h"
#include "EMGaussian.h"

using namespace std;

namespace Terran {

Cluster::Cluster(const vector<vector<double> > &data, const vector<double> &period, const vector<vector<Param> > &initialParams) : 
    dataset_(data),
    period_(period),
    paramset_(initialParams),
    partitions_(paramset_.size()) {

}

int Cluster::getNumDimensions() const {
    return paramset_.size();
}

int Cluster::getNumPoints() const {
    return dataset_.size();
}

bool Cluster::isPeriodic(int dimension) const {
    if(period_[dimension] == 0)
        return true;
    else
        return false;
}

double Cluster::getPeriod(int dimension) const {
    if(!isPeriodic(dimension)) {
        stringstream msg;
        msg << "getPeriod() exception, dimension " << dimension << " is not periodic!";
        throw(std::runtime_error(msg.str()));
    }
    return period_[dimension];
}

vector<double> Cluster::getPoint(int n) const {
    if(n >= getNumPoints()) {
        throw(std::runtime_error("Cluster::getPoint() out of bounds!"));   
    }
    return dataset_[n];
}

vector<double> Cluster::getDimension(int d) const {
    if(d >= getNumDimensions()) {
        throw(std::runtime_error("Cluster::getDimension() out of bounds!"));   
    }
    vector<double> data(getNumPoints());
    for(int i=0; i<getNumPoints(); i++) {
        data[i] = dataset_[i][d];
    }
    return data;
}

vector<Param> Cluster::getParameters(int dimension) const {
    return paramset_[dimension];
}

void Cluster::partition(int dimension, double threshold) {
    vector<Param> params = paramset_[dimension];
    vector<double> partition;
    if(isPeriodic(dimension)) {
        const double period = period_[dimension];
        MethodsPeriodicGaussian mpg(params, period);
        vector<double> minima = mpg.findMinima();
        for(int i=0; i < minima.size(); i++) {
            double val = periodicGaussianMixture(params, minima[i], period);
            if(val < threshold) {
                partition.push_back(minima[i]);
            }
        }
    } else {
        // implement code to partition non periodic mixture models

    }
    partitions_[dimension] = partition;
}

void Cluster::setPartitions(int dimension, const vector<double> &p) {
    if(dimension > getNumDimensions() - 1) {
        throw(std::runtime_error("Dimension out of bounds\n"));
    }
    partitions_[dimension] = p;
}

vector<int> Cluster::cluster() {

    // todo: parallelize on multiple threads
    /*
    for(int d = 0; d < getNumDimensions(); d++) {
        if(paramset_[d].size() == 0) {    
            stringstream ss;
            ss << "Fatal! Parameters not set for dimension " << d <<endl;
            throw(std::runtime_error(ss.str()));
        }
        // optimize the parameters
        optimizeParameters(d);
        // partition each dimension using parameters optimized by EM
        partition(d, 0.05);
    }
    */

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

void Cluster::optimizeParameters(int d) {
    vector<double> data;
    for(int n = 0; n < getNumPoints(); n++) {
        data.push_back(dataset_[n][d]);
    }
    // load initial set of parameters
    const vector<Param> &initialParams = paramset_[d];

    if(isPeriodic(d)) {
        if(paramset_[d].size() == 0) {
            //EMPeriodicGaussian epg(data, period_[d]);

            //epg.multiAdaptiveRun()
        } else {
            EMPeriodicGaussian epg(data, initialParams, period_[d]);
            epg.run();
        }
        //paramset_[d] = epg.getParams();
    } else {
        //EMGaussian eg(data, initialParams);
        if(paramset_[d].size() == 0) {
        
        } else { 
            //eg.run();
        }
        //paramset_[d] = eg.getParams();
    }
}

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
