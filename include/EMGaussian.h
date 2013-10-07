#ifndef EM_GAUSSIAN_H
#define EM_GAUSSIAN_H

#include "EM.h"

namespace Terran {

// Canonical Expectation Maximization of Gaussian Mixture Models
class EMGaussian : public EM {
public:
    EMGaussian(const std::vector<double> &data);
    EMGaussian(const std::vector<double> &data, const std::vector<Param> &params);
    ~EMGaussian();

    void MStep();

private:

    void mergeParams();

    double qkn(int k, int n) const;

    double domainLength() const;

};

}

#endif
