#ifndef PARTITIONER_H_
#define PARTITIONER_H_

#include <vector>

// a partition tool takes a 1D dataset and finds the "dips"

// currently this PartitionTool needs to be able to provide a function allowing for some
// form of visualization via an evaluate function. 

// all Partitioners must act on some marginalized data
namespace Terran {

class Partitioner {

public:

    Partitioner(const std::vector<double> dataset, double period) : dataset_(dataset), period_(period) {};
    virtual ~Partitioner() {};

    virtual std::vector<double> partition() = 0;
    
protected:

    const std::vector<double> dataset_;
    double period_;

};

}

#endif