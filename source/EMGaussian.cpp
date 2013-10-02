#include "EMGaussian.h"
#include "MathFunctions.h"

#include <limits>

namespace Terran {

EMGaussian::EMGaussian(const std::vector<double> &data) : 
    EM(data) {

}

EMGaussian::EMGaussian(const std::vector<double> &data, const std::vector<Param> &params) : 
    EM(data, params) {

}

EMGaussian::~EMGaussian() { 

}

void EMGaussian::MStep() {
    for(int k=0; k<params_.size(); k++) {
		//cout << k << endl;
        // Compute new mean
        double numeratorSum = 0;
        double denominatorSum = 0;
        for(int n=0; n<data_.size(); n++) {
            numeratorSum += pikn_[n][k]*data_[n];
        }
        for(int n=0; n<data_.size(); n++) {
            denominatorSum += pikn_[n][k];
        }
        params_[k].u = numeratorSum / denominatorSum;

        // Compute new standard deviation
        numeratorSum = 0;
        for(int n=0; n<data_.size(); n++) {
            double dx = data_[n]-params_[k].u;
            numeratorSum += pikn_[n][k]*dx*dx;
        }
        params_[k].s = sqrt(numeratorSum / denominatorSum);

        // Compute new probability
        params_[k].p = denominatorSum / data_.size(); 
    }
}

double EMGaussian::qkn(int k, int n) const {
   double pk = params_[k].p;
   double uk = params_[k].u;
   double sk = params_[k].s;    
   double xn = data_[n];
   return pk * gaussian(uk, sk, xn);
}

double EMGaussian::domainLength() const {
    double min =  numeric_limits<double>::max();
    double max = -numeric_limits<double>::max();
    for(int i=0; i < data_.size(); i++) {
        min = (data_[i] < min) ? data_[i] : min;
        max = (data_[i] > max) ? data_[i] : max;
    }
    return max-min;
}

}
