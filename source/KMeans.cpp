#include "KMeans.h"
#include <complex>
#include <limits>

using namespace std;

KMeans::KMeans(const std::vector<std::vector<double> > &dataset, const std::vector<bool> &isPeriodic) :
    dataset_(dataset), isPeriodic_(isPeriodic) {

}

// return a D sized vector given a set of point indices
vector<double> KMeans::computeMean(const vector<int> &points, const vector<bool> &isPeriodic) {

    vector<double> center;
    // for each dimension
    for(int d = 0; d < isPeriodic.size(); d++) {
        // for each point
        double mean = 0;
        if(isPeriodic[d]) {
            double real = 0;
            double imag = 0;
            double coord = dataset_[points[n]][d];
            for(int n = 0; n < points.size(); n++) {
                real += cos(coord);
                imag += sin(coord);
            }
            complex z(real, imag);
            z = z / points.size();
            mean = arg(z);
        } else {
            for(int n = 0; n < points.size(); n++) {
                mean += coord;
            }
            mean = mean / points.size();
        }
        center[d] = mean;
    }
    
    return center;  
}



// compute the distance between two points

// calculates the shortest distance between two points given a period
inline double periodicDifference(double x1, double x2, double period) {
    double diff = x1-x2;
    diff -= floor(diff/period+0.5)*period;
    return diff;
}

double KMeans::distance(const vector<double> &p1, const vector<double> &p2) {
    double sum = 0;
    if(p1.size() != p2.size()) {
        throw(std::runtime_error("KMeans::distance() - p1.size() != p2.size() "));
    }
    for(int d = 0; d < p1.size(); d++) {
        if(isPeriodic_[d]) {
            double diff = periodicDifference(p1[d], p2[d], 2*PI);
            sum += diff*diff;
        } else {
            double diff = p1[d]-p2[d];
            sum += diff*diff;
        }
    }
    return sqrt(sum);
}

// assign dataset into K groups based on current centers 
// note: can be sped up using triangle inequality presumably (periodic triangle inequality)?
vector<vector<int> > KMeans::assign(const vector<vector<double> > &centers) {
    
    vector<vector<int> > groups(centers.size());
    
    for(int n=0; n < dataset_.size(); n++) {
        double smallest_distance = numeric_limits<double>::max;
        double smallest_index;
        for(int k = 0; k < centers.size(); k++) {
            distance = distance(dataset_[n], centers[k]);
            if( distance < smallest_distance) {
                smallest_distance = distance;
                smallest_index = k;
            }
        }
        groups[k].push_back(n);
    }

    return groups;
    
}

// of size K x D
vector<int> KMeans::run(const vector<const vector<double> > &initialCenters) {

    vector<vector<double> > centers(initialCenters);
        
    
    for(int i=0; i < centers.size(); i++) {
        
    }

}