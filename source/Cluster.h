#ifndef TERRAN_H_
#define TERRAN_H_

#include <vector>
#include "Param.h"

/* Main interface class of the Cluster clustering algorithm */

// Usage:
// 1. Invoke EM on each degree of freedom to partition the domain into disjoint intervals.
// 2. Assign each point to a bucket using a map
// 3. Return a clustering based on the map
// 4. (Optional) Repeat 1-3 for each cluster until each degree of freedom has exactly 1 partition//    This is done by invoking the Cluster class recursively until no more DOFs can be found

// The Cluster class should not be responsible for doing the actual EM, this overlaps with the duties of the various EM classes.

// The resulting tree might look like:
//                     *
//                     0
//           __________|_________*
//           0         1         2
//                           *___|___*
//                           0   1   2
//                          / \     / \
//                         0   1   0   1
// 
// The leaves of tree are clusters

namespace Terran {

class Cluster {

public:

    Cluster(const std::vector<std::vector<double> > &data, const std::vector<double> &period, const std::vector<std::vector<Param> > &initialParams);

    ~Cluster();

    // returns number of dimensions in the dataset
    int getNumDimensions() const;

    // returns number of points in the dataset
    int getNumPoints() const;

    // returns true if dimension d is periodic
    bool isPeriodic(int d) const;

    // if dimension is periodic, returns the period of the dimension d
    // else an exception is thrown
    double getPeriod(int d) const;

    // return point n
    std::vector<double> getPoint(int n) const;

    // return marginalized values for dimension d
    std::vector<double> getDimension(int d) const;

    // returns the set of parameters in the mixture model fitting dimension d
    std::vector<Param> getParameters(int d) const;

    // set the parameters of the mixture model for dimension d
    void setParameters(int d, const std::vector<Param> &params);

    // set the partitions of dimension d
    void setPartitions(int d, const std::vector<double> &p);

    // returns an assignment of points into clusters
    // each dimension must have been partitioned either by means of:
    // setPartitions() or invoking partition()
    std::vector<int> cluster();

    // simple algorithm to optimize the parameters on a given dimension d using EM

    // parameters given, no multi run, no adaptive run [OK]
    // no parameters given, multi run, adaptive run [OK]

    // if no parameters are given for dimension, then multiAdaptive run is used
    // if parameters are explicitly given, no multiRun is done, though adaptiveRun is still
    // an option

    void optimizeParameters(int d);


private:

    // partition the dimension into disjoint intervals
    void partition(int dimension, double threshold);

    // assigns a given point to a bucket
    std::vector<short> assign(int point) const;

    // points are stored in dataset, size N x D
    // change to array of array later for performance
    const std::vector<std::vector<double> > dataset_;

    // describes the domain of each dimension
    // period of 0 indicates the dimension is not periodic
    const std::vector<double> period_;

    // the mixture model used to describe each DOF, size D x variable
    std::vector<std::vector<Param> > paramset_;

    // disjoint partitions of the domain
    std::vector<std::vector<double> > partitions_;
    
};

} // namespace Terran

#endif
