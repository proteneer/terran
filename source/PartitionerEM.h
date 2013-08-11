#ifndef PARTITIONER_EM_H_
#define PARTITIONER_EM_H_

#include "Partitioner.h"
#include "EM.h"

namespace Terran {

class PartitionerEM : public Partitioner {

public:
    
    PartitionerEM(const std::vector<double> &dataset, double period);
    ~PartitionerEM();
    // default method uses multiAdaptiveRun in EM.
    void optimizeParameters();

    double evaluate(double x) const;

    std::vector<double> partition();
    
    EM& getEM();

    bool isPeriodic() const;

private:
    
    std::vector<double> findLowMinima() const;

    // Multi adaptive run parameters
    // TODO: setters/getters for all these
    double MARThreshold_;
    int    MARNumParams_;
    int    MARNumTries_;

    // Minima whose value is less than partitionCutoff_ is 
    // considered to be a partition point
    double partitionCutoff_;

    double period_;
    EM* em_; 

};

}

#endif