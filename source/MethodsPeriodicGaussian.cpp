#include "MethodsPeriodicGaussian.h"
#include <stdexcept>
#include <assert.h>
#include <iostream>
#include <fstream>
#include <algorithm>  

using namespace std;

namespace Terran {

static void plotPeriodicGaussian(const vector<Param> &params, double period, int images) {
    ofstream mixture("mixture.dat");
    for(double xn = -period/2; xn < period/2; xn += 0.01) {
        mixture << xn << " " << periodicGaussianMixture(params, xn, period, images) << endl;
    }
    ofstream mixtureDx("mixtureDx.dat");
    for(double xn = -period/2; xn < period/2; xn += 0.01) {
        mixtureDx << xn << " " << periodicGaussianMixtureDx(params, xn, period, images) << endl;
    }
}

MethodsPeriodicGaussian::MethodsPeriodicGaussian(const vector<Param> &params, 
    double period) : Methods(params), period_(period) {

}

vector<double> MethodsPeriodicGaussian::findMaxima() const {
    vector<double> maxima;
    for(int k=0; k<params_.size(); k++) {
        double x_old = params_[k].u;
        bool found = false;
        const double delta = 1e-4;
        int iteration = 0;
        while(!found) {
            iteration++;
            if(iteration >= 1e3) {
                throw(std::runtime_error("MethodsPeriodicGaussian::findMaxima() - maximum iteration count reached!"));
            }
            double x_new = x_old - periodicGaussianMixtureDx(params_, x_old, period_)/periodicGaussianMixtureDx2(params_, x_old, period_);
            
            x_new = normalize(x_new);
            
            if(fabs(periodicGaussianMixtureDx(params_, x_new, period_)) < 1e-5)
                found = true;
            x_old = x_new;
        }
        
        // see if found point is a maxima
        double testLeft = normalize(x_old-delta);
        double testRight = normalize(x_old+delta);
        if(periodicGaussianMixtureDx(params_, testLeft, period_) > 0 && periodicGaussianMixtureDx(params_, testRight, period_) < 0) {
            bool skip = false;
            // see if found point has been found before
            for(int i=0; i<maxima.size(); i++)
                if(fabs(x_old-maxima[i]) < 1e-3)
                    skip = true;
            if(!skip) 
                maxima.push_back(x_old);
        }
    }
    return maxima;
}

vector<double> MethodsPeriodicGaussian::findMinima() const {
    vector<double> minima;
    vector<double> maxima = findMaxima();
    assert(maxima.size() > 0);
    sort(maxima.begin(), maxima.end());
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
            if(iteration > 1e5) {
                throw(std::runtime_error("Error: findMinima maximized number of iterations reached."));
            }
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
            my = periodicGaussianMixtureDx(params_, mk, period_);
            if(my < 0) 
                ak = mk;
            else
                bk = mk;
        } while (fabsp(ak,bk,period_) > 1e-8);
        minima.push_back(mk);
    }
    return minima;
}

}
