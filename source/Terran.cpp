#include "Terran.h"
#include "MethodsPeriodicGaussian.h"
#include "EMPeriodicGaussian.h"
#include "EMGaussian.h"

using namespace std;

Terran::Terran(const std::vector<const std::vector<double> > &data, const std::vector<double> &period, const std::vector<std::vector<Param> > &initialParams) : 
    data(data),
    period(period),
    params(initialParams) {

}

bool Terran::isPeriodic(int dimension) const {
    if(period[dimension] == 0)
        return true;
    else
        return false;
}

double Terran::getPeriod(int dimension) const {
    if(period[dimension] == 0) {
        throw(std::runtime_error("getPeriod() exception, dimension is not periodic!"));
    return period[dimension];
}

double Terran::getImages(int dimension) const {
    if(!isPeriodic(dimension))
        throw(std::runtime_error("getImages() exception, dimension is not periodic!"));
    return period[dimension];
}

void Terran::setImages(int dimension, int numImages) {
    if(!isPeriodic(dimension))
        throw(std::runtime_error("setImages() exception, dimension is not periodic!"));
    images[d] = numImages;
}

Terran::vector<double> partition(double threshold, int dimension) {
    vector<double> params = models_[dimension].params;
    vector<double> partition;
    try {
        const MixtureModelPeriodic &mmp = dynamic_cast<const MixtureModelPeriodic &>(models_[i]);
        double period = mmp.period;
        int images = mmp.images;
        MethodsPeriodicGaussian mpg(params, period, images);
        vector<double> minima = mpg.findMinima();
        for(int i=0; i<minima.size(); i++) {
            double val = periodicGaussianMixture(params, minima[i], period, images);
            if(val < threshold) {
                partition.push_back(minima[i]);
            }
        }
    } catch(const std::bad_cast &bc) {
        // implement code to partition non periodic mixture models

    }

    return partition;
}

void Terran::optimizeParameters(int d) {
    vector<double> data;
    for(int n = 0; n < getNumPoints(); n++) {
        data.push_back(dataset[n][d]);
    }
    
    // load initial set of parameters
    const vector<data> &initialParams = paramsset_[d];

    if(isPeriodic(d)) {
        em = new EMPeriodicGaussian epg(data, params, period[d], images[d]);
    } else {
        em = new EMGaussian eg(data, params);
        eg.run();
    }

    paramsset_[d] = EM.getParams();
}

vector<short> Terran::assign(int pointIndex) {
    assert(point.size() == intervals.size());
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

vector<int> Terran::cluster() const{

    // todo: parallelize on multiple threads
    vector<vector<double> > intervals;
    for(int d = 0; d < getNumDimensions(); d++) {
        // optimize the parameters
        optimizeParameters(d);

        // partition each dimension using parameters optimized by EM
        intervals_.push_back(partition(d));
    }

    // assign each point to a bucket
    map<vector<short>, vector<int> > clusters;
    for(int n = 0; n < getNumPoints(); n++) {
        vector<short> bucket = assign(n);
        clusters[bucket].push_back(n);
    }

    // loop over the clusters and get the points in each
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
