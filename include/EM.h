/****************************
 Copyright Yutong Zhao 2013
     Licensed under GPL3

     proteneer@gmail.com
****************************/

#ifndef EM_H_
#define EM_H_

#include <vector>
#include <stdexcept>
#include "Param.h"

// Abstract Expectation Maximization class for Gaussian-like mixture models that 
// optimize a set of initial parameters given a dataset. Concrete classes reimplement 
// the M step as needed. 
//
// Ref 1. Estimating Gaussian Mixture Densities with EM - A Tutorial, Carlo Tomasi
// 
// In general, the size of the data need not be exceedingly large provided the underlying
// data is well distributed. 

namespace Terran {

class EM { 

    public:

        EM(const std::vector<double> &data, const std::vector<Param> &params);
        EM(const std::vector<double> &data);
        virtual ~EM();

        // Set parameters
        void setParameters(const std::vector<Param> &input);

        // Get parameters
        std::vector<Param> getParams() const;

        // Get number of data points
        int getDataSize() const;

        // Compute the log likelihood given current parameters
        double getLikelihood() const;

        // Set maximum number of steps in any given EM run
        void setMaxSteps(int maxSteps);

        // Get the maximum number of steps
        int getMaxSteps() const;

        // Set the tolerance for the log likelihood convergence criteria
        void setTolerance(double tol);

        // Return the tolerance
        double getTolerance() const;

        // Runs the EM algorithm.
        // Notes:
        // -Requires parameters to be set explicitly
        // -Returns true if converged, false otherwise
        bool run();

        // Runs EM by guessing the number of initial parameters
        // Notes:
        // -Does not require parameters to be set explicitly
        // -Returns true if converged, false otherwise
        bool simpleRun(unsigned int numParams);

        // Compute the Expectation based on current parameters
        void EStep();

        // Maximize the Expectation by tuning parameters
        virtual void MStep() = 0;

    protected:

        // The sum over k for each n in p(k|n) should be 1
        void testIntegrity() const;
        
        const std::vector<double> &data_;
        std::vector<Param> params_;

        // pikn_ is a matrix of conditional probabilities: 
        // that given a point n was observed, it came from 
        // component k, ie. p(k|n) during iteration i
        // this is updated during the E-step
        std::vector<std::vector<double> > pikn_;

    private:

        // maximum number of steps in each EM run
        int maxSteps_;

        // tolerance threshold
        double tolerance_;

        // Used by the EStep to compute the expectation
        virtual double qkn(int k, int n) const = 0;

        // Estimate the domain size
        virtual double domainLength() const = 0;

        // Merge excessive parameters
        virtual void mergeParams() = 0;

        // this function cleans up weird parameters
        // 1. standard deviation too big/small
        // 2. probability weight too low (<0.01)
        bool cleanParameters();

};

}
#endif
