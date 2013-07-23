#include "Param.h"
#include <vector>
#include <map>
#include <assert.h>

using namespace std;

// intervals is of size D x # of buckets in the interval
vector<short> assignPointToBucket(const vector<double> point, const vector<bool> &isPeriodic, vector<vector<double> > &intervals) {

    assert(point.size() == intervals.size());

    vector<short> bucket(point.size());
    // for each dimension in the point
    for(int i=0; i<point.size(); i++) {
        // for each set of intervals
        for(int j=0; j<intervals[i].size(); j++) {
            if(point[i] < intervals[i][j]) {
                bucket[i] = j;
                break;
            }

            if(isPeriodic[i]) {
                bucket[i] = 0;
            } else {
                bucket[i] = j+1;
            }
        }
    }
    return bucket;
}

// data - N x D collection of points
// isPeriodic - D x 1 vector denoting if dimension d is periodic
// intervals - D x # of buckets
vector<int> assignDataset(const vector<vector<double> > &data, vector<bool> isPeriodic, vector<vector<double > > &intervals) {

    // maps a bucket ID to a set of points
    map<vector<short>, vector<int> > bucketCount;

    for(int i=0; i < data.size(); i++) {
        const vector<double> &point = data[i];
        bucketCount[assignPointToBucket(point, isPeriodic, intervals)].push_back(i);
    }

    int clusterIndex = 0;
    vector<int> cluster(data.size());
    for(map<vector<short>, vector<int> >::const_iterator it = bucketCount.begin(); it != bucketCount.end(); it++) {
        vector<int> points = it->second;
        for(int i=0; i<points.size(); i++) {
            cluster[points[i]] = clusterIndex;
        }
        clusterIndex++;
    }

}
