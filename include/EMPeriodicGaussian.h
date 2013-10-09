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

		void EStep();

        void MStep();

    private:

		void initializePink();

		// pinky_[n][k][r] - where r can go from -6 to 6 -> mapped to 0 to 12
		std::vector<std::vector<std::vector<double> > > pinkr_;
        
        void mergeParams();

        double domainLength() const;

        double qkn(int k, int n) const; 

        // Simplified versions of the lambda derivative in order to find roots
		/*
        double dlds(double sk, int k) const;
        double dldu(double uk, int k) const;
		*/

        double period_;
};

}
#endif
