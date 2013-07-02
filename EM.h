/****************************
 Copyright Yutong Zhao 2013
     Licensed under GPL3

     proteneer@gmail.com
****************************/

#ifndef EM_H_
#define EM_H_

#include <vector>
#include "Param.h"

// Abstract Expectation Maximization class for Gaussian-like mixture models that 
// optimize a set of initial parameters given a dataset. Concrete classes reimplement 
// the E and M step as needed.
class EM { 
    public:
        virtual ~EM();

        // Run the EM algorithm
    	void run(int numSteps);

        // Set parameters
        void setParams(std::vector<Param> &input);

        // Get parameters
        std::vector<Param> getParams() const;

        int getDataSize() const {
            return data_.size();
        }

        // EM routines
        virtual void EStep() = 0;
        virtual void MStep() = 0;

    private:
        // The sum over k for each n in p(k|n) should be 1
        // 
        void testIntegrity() const;

        // An extension, as defined in the paper
        double qkn(int k, int n) const;

        std::vector<Param> params_;
        const std::vector<double> data_;

        // pikn_ is a matrix of conditional probabilities: 
        // that given a point n was observed, it came from 
        // component k, ie. p(k|n) during iteration i
        // this is updated during the E-step
        std::vector<std::vector<double> > pikn_;
};
#endif
