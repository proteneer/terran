#ifndef EM_GAUSSIAN_H
#define EM_GAUSSIAN_H

#include "EM.h"

namespace Terran {

// Canonical Expectation Maximization of Gaussian Mixture Models
class TERRAN_EXPORT EMGaussian : public EM {
public:
    EMGaussian(const std::vector<double> &data);
    EMGaussian(const std::vector<double> &data, const std::vector<Param> &params);
    ~EMGaussian();

    void MStep();

	void EStep();
	
private:

	// initialize internals based on current data_, params_
	void initializePink();

	void destroyPink();

	// pikn_ is a matrix of conditional probabilities: 
    // that given a point n was observed, it came from 
    // component k, ie. p(k|n) during iteration i
    // this is updated during the E-step
    std::vector<std::vector<double> > pink_;

    void mergeParams();

    double qkn(int k, int n) const;

    double domainLength() const;

};

}

#endif
