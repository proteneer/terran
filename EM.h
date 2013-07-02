#ifndef _EM_HH_
#define _EM_HH_

#include <vector>
#include "Params.hh"

#define CALL_MEMBER_FN(object, ptm) ((object).*(ptm))

using namespace std;

struct Param {
    // probability weight of this mixture
    double p; 
    // mean
    double u;
    // standard deviation
    double s;
}

// Abstract Expectation Maximization class that optimize a set of initial parameters 
// given a dataset. Concrete classes reimplement the E and M step
class EM { 
    public:
        virtual ~EM();

        // Run the EM algorithm
    	void run(void numSteps);

        // Set parameters
        void setParams(vector<Param> &input);

        // Get parameters
        vector<Param> getParams() const;

        int getDataSize() const {
            return data_.size();
        }

        // EM Routines
        virtual void EStep() = 0;
        virtual void MStep() = 0;

    private:
        // The sum over k for each n in p(k|n) should be 1
        void assert_pikn_integrity() const;

        // An extension, as defined in the paper
        double q_k_n(int k, int n) const;

        vector<Param> params_;
        const vector<double> data_;

        // p_i_k_n_ is a matrix of conditional probabilities: 
        // that given a point n was observed, it came from 
        // component k, ie. p(k|n) during iteration i
        // this is updated during the E-step
        vector<vector<double> > pikn_;

};
#endif
