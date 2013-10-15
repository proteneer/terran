#ifndef TERRAN_H_
#define TERRAN_H_

#include <vector>
#include <string>

#include "export.h"
#include "Param.h"
#include "Partitioner.h"

namespace Terran {
    
class TERRAN_EXPORT Cluster {

public:

	Cluster(const std::vector<std::vector<double> > &data, const std::vector<int> &period, Partitioner* partitioner); 
    
	Cluster(const std::vector<std::vector<double> > &data, const std::vector<int> &period);

    // todo: delete partitions!
    ~Cluster();

    // returns number of dimensions in the dataset
    int getNumDimensions() const;

    // returns number of points in the dataset
    int getNumPoints() const;

    // returns true if dimension d is periodic
    bool isPeriodic(int d) const;

    // return point n of length D
    std::vector<double> getPoint(int n) const;

    // return marginalized values for dimension d of length N
    std::vector<double> getDimension(int d) const;
    
    // returns the partition for dimension d
    std::vector<double> getPartition(int d) const;

    // set the partition info explicitly for dimension d
    void setPartition(int d, const std::vector<double> &p);

    // set the partition method. Currently supported: "EM", possible "KDE" in the future
    void setPartitionMethod(std::string method = "EM");

    // partition dimension d
    void partition(int d);

	// partition every dimension
	void partitionAll();

    // returns an assignment of points into clusters
    // each dimension must have been partitioned either by means of:
    // setPartitions() or invoking partition()
    std::vector<int> assign();

    void setSubsampleCount(int count);

    int getSubsampleCount() const;

private:
	
	void initialize();

    // assigns a given point to a bucket
    std::vector<short> findBucket(int point) const;

    std::string partitionMethod_;

    // number of points used to subsample
    int subsampleCount_;

    // points are stored in dataset, size N x D
    // change to array of array later for performance
    const std::vector<std::vector<double> > dataset_;

	// of size d, partitionFlag_[d] is true if dimension d
	// has been partitioned, false otherwise.
	// all dimensions must have been partitioned at least once
	std::vector<bool> partitionFlag_;

    // describes the domain of each dimension
    // period of 0 indicates the dimension is not periodic

	// enforce period of 2PI.
    const std::vector<int> period_;

    // disjoint partitions of each domain
    std::vector<std::vector<double> > partitions_;  

	Partitioner* partitioner_;

};

} // namespace Terran

#endif
