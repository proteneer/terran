#include "PartitionPeriodicGaussian.h"
#include <stdexcept>
#include <assert.h>

using namespace std;

PartitionPeriodicGaussian::PartitionPeriodicGaussian(const vector<Param> &params, 
    double period, int images) : Partition(params), period_(period), images_(images) {

}

void PartitionPeriodicGaussian::partition(double threshold) const {

}

vector<double> PartitionPeriodicGaussian::findMaxima() const {
    vector<double> maximas;
    for(int k=0; k<params_.size(); k++) {
        double xn_old = params_[k].u;
        bool found = false;
        const double delta = 1e-4;

        // careful: change this to a long if need to be >2e9
        int iteration = 0;
        while(!found) {
            iteration++;
            // Saddle points suck. If the simple algorithm doesn't converge
            // it's probably due to a saddle point. 
            if(iteration >= 1e7) {
                throw(std::runtime_error("findPeriodicMaxima: maximum iteration count reached!"));
            }
            double xn_new = xn_old + delta*periodicGaussianMixtureDx(params_, xn_old, period_, images_);
            xn_new = normalize(xn_new);
            if(fabsp(xn_new, xn_old, period_) < 1e-8)
                found = true;
            xn_old = xn_new;
        }
        bool skip = false;
        // if two gaussians are too close to each other, discard it
        for(int i=0; i<maximas.size(); i++)
            if(fabsp(xn_old, maximas[i], period_) < 1e-3)
                skip = true;
        if(!skip) 
            maximas.push_back(xn_old);
    }
    return maximas;

}

vector<double> PartitionPeriodicGaussian::findMinima() const {
    vector<double> minima;
    vector<double> maxima = findMaxima();
    assert(maxima.size() > 0);
    sort(maxima.begin(), maxima.end());
    const double tol = 1e-6;
    const double perturb = 1e-2;
    for(int k=0; k<maxima.size(); k++) {
        double ak = maxima[k];
        double bk = maxima[(k+1)%maxima.size()];
        ak = normalize(ak+perturb);
        bk = normalize(bk-perturb);
        double mk = 0;
        double my = 0;
        int iteration = 0;
        do {
            if(iteration > 1e6)
                throw(std::runtime_error("Error: findPeriodicMinimumBS maximized number of iterations reached."));
            else
                iteration++;
            
            //suppose we have the maximas:
            // a b c d, where a<b<c<d 
            // there are minimas in (a,b) (b,c) (c,d) (d,a)
            // the if statement takes care of the (d,a) case
            if(bk < ak)
                mk = normalize((ak+period_+bk)/2.0);
            else
                mk = normalize((ak+bk)/2.0);
            my = periodicGaussianMixtureDx(params_, mk, period_, images_);
            if(my < 0) 
                ak = mk;
            else
                bk = mk;
        //} while (fabs(my) > 1e-9);
        } while (fabsp(ak,bk,period_) > 1e-8);
        minima.push_back(mk);
    }
    return minima;
}
