#ifndef EM_GAUSSIAN_H
#define EM_GAUSSIAN_H

#include "EM.h"

namespace Terran {

// Canonical Expectation Maximization of Gaussian Mixture Models
class EMGaussian : public EM {
public:
    EMGaussian(const std::vector<double> &data, const std::vector<Param> &params);
    ~EMGaussian();

    void MStep();

private:

    double qkn(int k, int n) const;

    // given min and max of the domain
    // samples a set of points in the interval [min, max]
    std::vector<double> sampleDomain(int count) const;

    // samples a set of poitns in the interval [left, right]
    std::vector<double> sampleDomain(int count, double left, double right) const;

};

}

#endif
