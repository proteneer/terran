#ifndef EM_PERIODIC_GAUSSIAN_H
#define EM_PERIODIC_GAUSSIAN_H

#include "EM.h"
#include "MathFunctions.h"

namespace Terran {

// Expectation Maximization of Periodic Gaussian Mixture Models
class TERRAN_EXPORT EMPeriodicGaussian : public EM {
    public:
        
        explicit EMPeriodicGaussian(const std::vector<double> &data, const std::vector<Param> &params, double period);
        
        explicit EMPeriodicGaussian(const std::vector<double> &data, double period);
        ~EMPeriodicGaussian();

		void EStep();

        void MStep();

    private:

		void initializePink();

		void destroyPink();

		std::vector<std::vector<std::vector<double> > > pinkr_;
        
        void mergeParams();

        double domainLength() const;

        double qkn(int k, int n) const; 

        double period_;
};

}
#endif
