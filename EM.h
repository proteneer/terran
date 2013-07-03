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
#include "MathFunctions.h"

// Abstract Expectation Maximization class for Gaussian-like mixture models that 
// optimize a set of initial parameters given a dataset. Concrete classes reimplement 
// the E and M step as needed. 
//
// Ref 1. Estimating Gaussian Mixture Densities with EM - A Tutorial, Carlo Tomasi
class EM { 
    public:
        EM(const std::vector<double> &data, const std::vector<Param> &params) : 
            data_(data),
            params_(params),
            pikn_(data.size(), std::vector<double>(params.size(),0)) {
        
        }

        virtual ~EM();

        // Set parameters
        void setParams(const std::vector<Param> &input) {
            params_ = input;
        }

        // Get parameters
        std::vector<Param> getParams() const {
            return params_;
        }

        // Get number of data points
        int getDataSize() const {
            return data_.size();
        }

        // Run the EM algorithm
    	virtual void run(int maxSteps, double tolerance) = 0;

        // EM routines
        void EStep() {
            for(int n=0; n<data_.size(); n++) {
                double sum = 0;
                for(int k=0; k<params_.size(); k++) {
                    sum += qkn(k,n);
                }
                for(int k=0; k<params_.size(); k++) {
                    pikn_[n][k] = qkn(k,n)/sum;
                }
            }
            testIntegrity();
        }

        virtual void MStep() = 0;

    private:

        virtual double qkn(int k, int n) const = 0;

    protected:
        // The sum over k for each n in p(k|n) should be 1
        void testIntegrity() const {
            double sum = 0;
            for(int k=0; k<params_.size(); k++) {
                for(int n=0; n<data_.size(); n++) {
                    sum += pikn_[n][k];   
                }
                
                if(fabs(sum-1.0) < 0.01) {
                    throw(std::runtime_error("pikn no longer sums to 1"));
                }
            }
        };

protected:
        const std::vector<double> data_;
        std::vector<Param> params_;

        // pikn_ is a matrix of conditional probabilities: 
        // that given a point n was observed, it came from 
        // component k, ie. p(k|n) during iteration i
        // this is updated during the E-step
        std::vector<std::vector<double> > pikn_;
};
#endif
