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

    // default method uses multiAdaptiveRun in EM.
    void optimizeParameters();

    // invokes optimizeParameters and findLowMinima
    std::vector<double> partition();
    
    EM& getEM();
    
private:

    bool isPeriodic() const;
    
    std::vector<double> findLowMinima() const;

    // Multi adaptive run parameters
    // TODO: setters/getters for all these
    double MARThreshold_;
    int    MARNumParams_;
    int    MARNumTries_;

    // Minima whose value is less than partitionCutoff_ is 
    // considered to be a partition point
    double partitionCutoff_;
    EM* em_; 

};

}

#endif