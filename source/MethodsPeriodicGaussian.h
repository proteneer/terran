#include "MathFunctions.h"
#include "Methods.h"

namespace Terran {

class MethodsPeriodicGaussian : public Methods {

public:

    explicit MethodsPeriodicGaussian(const std::vector<Param> &params,
        double period, int images);

    // Partition the domain into disjoint intervals
    std::vector<double> partition(double threshold) const;

    // Find the maxima for a periodic gaussian mixture model
    std::vector<double> findMaxima() const;

    // Find the minima for a periodic gaussian mixture model
    std::vector<double> findMinima() const;

private:

    const double period_;
    const int images_;
};

}
