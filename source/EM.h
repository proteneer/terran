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

        // Run the canonical EM algorithm 
        // Returns true if executed successfully, false otherwise
        // bool run(int maxSteps = 100, double tolerance = 0.1, const vector<Param> &params_);
        // Requires parameters to be set explicitly
    	bool run();

        // Run the adaptive EM algorithm
        // The number of parameters that overfit will be pruned away automatically
        // This is determined by the probability weight component p
        // Requires parameters to be set explicitly
        bool adaptiveRun(double cutoff);

        // Multiadaptive run generates a set of parameters by random repeatedly
        // and sees which results in the highest likelihood

        // make it possible to set the left and right of the sample domain
        void multiAdaptiveRun(double cutoff, int numParams, int numTries);

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

        // Generate random samples in the domain
        virtual std::vector<double> sampleDomain(int count) const = 0;

        // Generates random samples in the region [left, right]
        virtual std::vector<double> sampleDomain(int count, double left, double right) const = 0;
};

}
#endif
