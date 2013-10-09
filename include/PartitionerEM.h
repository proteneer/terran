#ifndef PARTITIONER_EM_H_
#define PARTITIONER_EM_H_

#include "Partitioner.h"
#include "EM.h"

namespace Terran {

/* Implements the Partitioner class using the expectation maximization algorithm.

   By default, EM::multiAdaptiveRun() is used to fit a set of Gaussian Mixture Models
   to the underlying data. The resulting minima of said curve is found, and the subset
   of minima whose probability is less than the cutoff is set as the partition.
   
   For more fined tuned control of how the EM should be done,
   such as setting initial gaussian mixture parameters, a reference to the underlying EM
   object can be obtained.
*/
class PartitionerEM : public Partitioner {

public:
    
    PartitionerEM(const std::vector<double> &dataset, double period);
    ~PartitionerEM();

    // Executes EM::multiAdaptiveRun()
    void optimizeParameters();

    // invokes optimizeParameters and findLowMinima
    std::vector<double> partition();
    
    EM& getEM();

    void setPartitionCutoff(double val);

    double getPartitionCutoff() const;

private:

    bool isPeriodic() const;
    
    std::vector<double> findLowMinima() const;

    // To do: establish left and right bounds for sampling the mean

    // During MAR, when a given parameter for some mixture model
    // has a weight < MARThreshold_, the parameter is purged
    double MARThreshold_;

    // The number of parameters to start with for MAR
    int    MARNumParams_;

    // Number of attempts
    int    MARNumTries_;

    // Minima whose value is less than partitionCutoff_ is 
    // considered to be a partition point
    double partitionCutoff_;
    EM* em_; 

};

}

#endif