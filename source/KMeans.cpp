#include "KMeans.h"
#include "MathFunctions.h"
#include <complex>
#include <limits>

using namespace std;

namespace Terran {

KMeans::KMeans(const std::vector<std::vector<double> > &dataset, const std::vector<bool> &isPeriodic) :
    dataset_(dataset), isPeriodic_(isPeriodic) {
    if(dataset_.size() == 0) {
        throw(std::runtime_error("KMeans::KMeans() - Empty dataset!"));
    }
    for(int n = 0; n < dataset_.size(); n++) {
        if(dataset_[n].size() != isPeriodic.size()) {
            throw(std::runtime_error("KMeans::KMeans() - dataset size does not match isPeriodic size"));
        }
    }
}

KMeans::~KMeans() {}

// return a D sized mean vector given a list of point indices
vector<double> KMeans::computeMean(const vector<int> &points) const {

    vector<double> center(isPeriodic_.size());
    // for each dimension
    for(int d = 0; d < isPeriodic_.size(); d++) {
        double mean = 0;
        // if the dimension is periodic we use directional statistics
        // to compute the periodic mean
        if(isPeriodic_[d]) {
            double real = 0;
            double imag = 0;
            for(int n = 0; n < points.size(); n++) {
                double coord = dataset_[points[n]][d];
                real += cos(coord);
                imag += sin(coord);
            }
            complex<double> z(real/points.size(), imag/points.size());
            mean = arg(z);
        } else {
            for(int n = 0; n < points.size(); n++) {
                double coord = dataset_[points[n]][d];
                mean += coord;
            }
            mean = mean / points.size();
        }
        center[d] = mean;
    }
    return center;  
}

double KMeans::distance(const vector<double> &p1, const vector<double> &p2) const {
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
 vector<int> KMeans::assignPointsToCenters(const vector<vector<double> > &centers) const {
    // partition N points into K groups of variable length
    vector<int> assignment(dataset_.size());
    for(int n=0; n < dataset_.size(); n++) {
        double smallest_distance = numeric_limits<double>::max();
        int smallest_index;
        for(int k = 0; k < centers.size(); k++) {
            double dist = distance(dataset_[n], centers[k]);
            if( dist < smallest_distance) {
                smallest_distance = dist;
                smallest_index = k;
            }
        }
        assignment[n] = smallest_index;
    }
    return assignment;
}


// of size K x D
std::vector<int> KMeans::run(const std::vector<std::vector<double> > &initialCenters) const {

    bool converged = false;  
    vector<int> assignment = assignPointsToCenters(initialCenters);

	int steps = 0;

    while(!converged && steps < 1e5) {
		// partition N points into K groups of variable size
        vector<vector<int> > groups(initialCenters.size()); 
        for(int n=0; n < assignment.size(); n++) {
            groups[assignment[n]].push_back(n);
        }
    
        // compute centroid of each group of points
        vector<vector<double> > centers(initialCenters);
        for(int k = 0; k < groups.size(); k++) {
            centers[k] = computeMean(groups[k]);
        }

        // compute new assignment based on these k groups
        vector<int> newAssignment = assignPointsToCenters(centers);

        if(newAssignment != assignment) {
            assignment = newAssignment;
        } else {
            converged = true;
        }
		steps++;
    }

	return assignment;
}

}