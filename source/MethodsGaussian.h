#ifndef METHODS_GAUSSIAN_H_
#define METHODS_GAUSSIAN_H_

#include "MathFunctions.h"
#include "Methods.h"

namespace Terran {

class MethodsGaussian : public Methods {

public:

    explicit MethodsGaussian(const std::vector<Param> &params);

    // Partition the domain into disjoint intervals
    // std::vector<double> partition(double threshold) const;

    // Find the maxima for a gaussian mixture model
    std::vector<double> findMaxima() const;

    // Find the minima for a gaussian mixture model
    std::vector<double> findMinima() const;

};

}

#endif