#ifndef EM_PERIODIC_GAUSSIAN_H
#define EM_PERIODIC_GAUSSIAN_H

#include "EM.h"
#include "MathFunctions.h"

namespace Terran {

// Expectation Maximization of Periodic Gaussian Mixture Models
class EMPeriodicGaussian : public EM {
    public:
        
        explicit EMPeriodicGaussian(const std::vector<double> &data, const std::vector<Param> &params, double period);
        
        explicit EMPeriodicGaussian(const std::vector<double> &data, double period);
        ~EMPeriodicGaussian();

        void MStep();

    private:
        // sample points in the periodic domain
        virtual vector<double> sampleDomain(int count) const;

        double qkn(int k, int n) const; 

        // Simplified versions of the lambda derivative in order to find roots
        double dlds(double sk, int k) const;
        double dldu(double uk, int k) const;

        double period_;
};

}
#endif
