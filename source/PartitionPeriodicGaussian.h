#include "Partition.h"
#include "MathFunctions.h"

class PartitionPeriodicGaussian : public Partition {

public:

    PartitionPeriodicGaussian(const std::vector<Param> &params,
        double period = 2*PI, int images = 10);

    // Partition the domain into disjoint intervals
    void partition(double threshold) const;

    // Find the maxima
    std::vector<double> findMaxima() const;

    // Find the minima
    std::vector<double> findMinima() const;

private:

    double period_;
    int images_;
};
