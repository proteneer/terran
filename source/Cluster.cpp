#include <stdexcept>
#include <sstream>

#include "Cluster.h"
#include "MethodsPeriodicGaussian.h"
#include "EMPeriodicGaussian.h"
#include "EMGaussian.h"

using namespace std;

namespace Terran {

Cluster::Cluster(const std::vector<const std::vector<double> > &data, const std::vector<double> &period, const std::vector<std::vector<Param> > &initialParams) : 
    dataset_(data),
    period_(period),
    paramset_(initialParams),
    images_(vector<int>(paramset_.size(),15)),
    partitions_(vector<vector<double> >(paramset_.size(),vector<double>())) {

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

double Cluster::getImages(int dimension) const {
    if(!isPeriodic(dimension)) {
        stringstream msg;
        msg << "getImages() exception, dimension " << dimension << " is not periodic!";
        throw(std::runtime_error(msg.str()));
    }
    return period_[dimension];
}

double Cluster::getParameters(int dimension) const {
    return paramset_[dimension];
}

void Cluster::setImages(int dimension, int numImages) {
    if(!isPeriodic(dimension)) {
        stringstream msg;
        msg << "setImages() exception, dimension " << dimension << " is not periodic!";
        throw(std::runtime_error(msg.str()));
    }
    images_[dimension] = numImages;
}

Cluster::vector<double> partition(int dimension, double threshold) {
    vector<double> params = paramset_[dimension];
    vector<double> partition;
    if(isPeriodic(d)) {
        const double period = period_[d];
        const int images = images_[d];
        MethodsPeriodicGaussian mpg(params, period, images);
        vector<double> minima = mpg.findMinima();
        for(int i=0; i < minima.size(); i++) {
            double val = periodicGaussianMixture(params, minima[i], period, images);
            if(val < threshold) {
                partition.push_back(minima[i]);
            }
        }
    } else {
        // implement code to partition non periodic mixture models

    }
    return partition;
}

vector<int> Cluster::run() const{

    // todo: parallelize on multiple threads
    vector<vector<double> > intervals;
    for(int d = 0; d < getNumDimensions(); d++) {
        // optimize the parameters
        optimizeParameters(d);
        // partition each dimension using parameters optimized by EM
        partitions_.push_back(partition(d));
    }

    // assign each point to a bucket
    // cannot parallelize easily
    map<vector<short>, vector<int> > clusters;
    for(int n = 0; n < getNumPoints(); n++) {
        vector<short> bucket = assign(n);
        clusters[bucket].push_back(n);
    }

    // loop over the buckets and get the points in each
    int clusterIndex = 0;
    vector<int> assignment(getNumPoints(),-1);
    for(map<vector<short>, vector<int> >::const_iterator it = clusters.begin();
        it != clusters.end(); it++) {
        vector<int> pointIndices = it->second;
        for(int j=0; j<pointIndices.size(); j++) {
            assignment[pointsIndices[j]] = clusterIndex;
        }
        clusterIndex++;
    }

    return assignment;
}

void Cluster::optimizeParameters(int d) {
    vector<double> data;
    for(int n = 0; n < getNumPoints(); n++) {
        data.push_back(dataset[n][d]);
    }
    // load initial set of parameters
    const vector<data> &initialParams = paramsset_[d];

    if(isPeriodic(d)) {
        EMPeriodicGaussian epg(data, initialParams, period[d], images[d]);
        epg.run();
        paramsset_[d] = epg.getParams();
    } else {
        EMGaussian eg(data, initialParams);
        eg.run();
        paramsset_[d] = eg.getParams();
    }
}

vector<short> Cluster::assign(int pointIndex) {
    vector<short> bucket(getNumDimensions());
    vector<double> point = dataset[pointIndex];

    // for each dimension in the point
    for(int d=0; d<getNumDimensions(); d++) {
        // loop over the intervals
        for(int j=0; j<intervals[d].size(); j++) {
            if(point[d] < intervals[d][j]) {
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
