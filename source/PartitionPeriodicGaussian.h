#include "Partition.h"
#include "MathFunctions.h"

class PartitionPeriodicGaussian : public Partition {
public:

    explicit PartitionPeriodicGaussian(const std::vector<Param> &params,
        double period = 2*PI, int images = 10);

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
