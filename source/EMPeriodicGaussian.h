#ifdef EM_PERIODIC_GAUSSIAN_H
#define EM_PERIODIC_GAUSSIAN_H

#include "EM.h"

// Expectation Maximization of Periodic Gaussian Mixture Models
class EMPeriodicGaussian : public EM {
    public:
        EMPeriodicGaussian(const std::vector<double> &data, const std::vector<Param> &params);
        ~EMPeriodicGaussian();

        void MStep();
    private:
        double qkn(int k, int n) const; 

        double dlds(double sk, int k) const;
        double dldu(double uk, int k) const;

        double period_;
        double images_;
}
