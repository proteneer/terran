#ifndef TERRAN_H_
#define TERRAN_H_

#include <vector>
#include <string>
#include "Param.h"
#include "Partitioner.h"
#include "PartitionerEM.h"

namespace Terran {
    
class Cluster {

public:

    Cluster(const std::vector<std::vector<double> > &data, const std::vector<double> &period);

    // todo: delete partitions!
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

    // return point n of length D
    std::vector<double> getPoint(int n) const;

    // return marginalized values for dimension d of length N
    std::vector<double> getDimension(int d) const;
    
    // returns the partition for dimension d
    std::vector<double> getPartition(int d) const;

    // set the partition info explicitly for dimension d
    void setPartition(int d, const std::vector<double> &p);

    void setPartitionMethod(std::string method = "EM");

    // partition dimension d
    void partition(int d);

    // returns an assignment of points into clusters
    // each dimension must have been partitioned either by means of:
    // setPartitions() or invoking partition()
    std::vector<int> cluster();

    void setSubsampleCount(int count);

    int getSubsampleCount() const;

private:

    // assigns a given point to a bucket
    std::vector<short> assign(int point) const;

    std::string partitionMethod_;

    // number of points used to subsample
    double subsampleCount_;

    // points are stored in dataset, size N x D
    // change to array of array later for performance
    const std::vector<std::vector<double> > dataset_;

    // describes the domain of each dimension
    // period of 0 indicates the dimension is not periodic
    const std::vector<double> period_;

    // disjoint partitions of each domain
    std::vector<std::vector<double> > partitions_;  


};

} // namespace Terran

#endif
