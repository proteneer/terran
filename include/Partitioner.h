#ifndef PARTITIONER_H_
#define PARTITIONER_H_

#include <vector>
#include "export.h"

/* Abstract Class

   Partitioner take a given set of 1d data, and a flag denoting the period of the data.
   If the period is zero, then the data is assumed to be non-periodic.

   Implementations of Partitioner must implement their own partition() functions. It is
   safe to make references to dataset_ and period_ as needed. Semantically, a partitioner
   can makes low density cuts in the dataset.
*/
namespace Terran {

class TERRAN_EXPORT Partitioner {

public:

    Partitioner(const std::vector<double> dataset, double period) : dataset_(dataset), period_(period) {};
    virtual ~Partitioner() {};

    virtual std::vector<double> partition() = 0;
    
protected:

    // The returned vector is defined as follows:
    // - the resulting partitions must be sorted
    // - empty vector implies no partitions
    // - otherwise, for a periodic system, the returned vector of size
    //   n is a partition into n parts, _0_|__1__|_0_
    //   for a nonperiodic system, the returned vector of size n
    //   is a partition into n+1 parts, _0_|__1__|_2_
    //   where | denote the partition element
    const std::vector<double> dataset_;
    double period_;

};

}

#endif