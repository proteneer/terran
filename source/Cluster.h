#ifndef TERRAN_H_
#define TERRAN_H_

#include <vector>
#include "MixtureModel.h"

/* Main interface class of the Cluster clustering algorithm */


// Usage:
// 1. Invoke EM on each degree of freedom to partition the domain into disjoint intervals.
// 2. Assign each point to a bucket using a map
// 3. Return a clustering based on the map
// 4. (Optional) Repeat 1-3 for each cluster until each degree of freedom has exactly 1 partition//    This is done by invoking the Cluster class recursively until no more DOFs can be found

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

    Cluster(const std::vector<const std::vector<double> > &data, const std::vector<double> &period, const std::vector<std::vector<Param> > &initialParams);

    ~Cluster();

    // returns number of dimensions in the dataset
    int getNumDimensions() const;

    // returns number of points in the dataset
    int getNumPoints() const;

    // returns true if dimension X is periodic
    bool isPeriodic(int dimension) const;

    // if dimension is periodic, returns the period of the dimension
    // else an exception is thrown
    double getPeriod(int dimension) const;

    // if dimension is periodic, returns the number of images of the dimension
    // else an exception is thrown
    double getImages(int dimension) const;

    // returns the set of parameters in the mixture model fitting dimension
    std::vector<Param> getParameters(int dimension) const;

    // set the parameters needed to approximate dimension
    void setParameters(int dimension, std::vector<Param> params);

    // set number of images to use
    void setImages(int dimension, int numImages);

    // returns an assignment of points into clusters
    std::vector<int> run();

    // run EM to optimize the parameters on a given dimension
    void optimizeParameters(int dimension);

private:

    // partition the dimension into disjoint intervals
    std::vector<double> partition(int dimension, double threshold);

    // assigns a given point to a bucket
    std::vector<short> assign(int point) const;

    // points are stored in dataset, size N x D
    const std::vector<std::vector<double> > dataset_;

    // describes the domain of each dimension
    // period of 0 indicates the dimension is not periodic
    const std::vector<double> period_;

    // the mixture model used to describe each DOF, size D x variable
    std::vector<std::vector<Param> > paramset_;

    // describes the number of images needed to approximate each dimension
    // this only needs to be mucked with if the standard deviation is large
    std::vector<int> images_;

    // disjoint partitions of the domain
    std::vector<std::vector<int> > partitions_;
    
}

} // namespace Terran

#endif
