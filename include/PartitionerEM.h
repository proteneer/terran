#ifndef PARTITIONER_EM_H_
#define PARTITIONER_EM_H_

#include "export.h"
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
class TERRAN_EXPORT PartitionerEM : public Partitioner {

public:
    
	PartitionerEM();
    
	~PartitionerEM();

	// From base-class
	void setDataAndPeriod(const std::vector<double> &data, bool isPeriodic);

    // invokes optimizeParameters and findLowMinima
    std::vector<double> partition();
    
	// get a curve of the model approximating the data in the range [left, right]
	void evaluateModel(std::vector<double> &x, std::vector<double> &y, int nsamples);

    EM& getEM();

    void setPartitionCutoff(double val);

    double getPartitionCutoff() const;

	void setInitialK(int count);

	int getInitialK() const;

private:
	// Executes EM::simpleRun()
    void optimizeParameters();

    std::vector<double> findLowMinima() const;

	// if this dataset is periodic
	bool isPeriodic_;

	// Num of components to use in simple run
	int initialK_;

    // Minima whose value is less than partitionCutoff_ is 
    // considered to be a partition point
    double partitionCutoff_;
    EM* em_; 

};

}

#endif