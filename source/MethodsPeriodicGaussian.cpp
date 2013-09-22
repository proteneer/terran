#include "MethodsPeriodicGaussian.h"
#include <stdexcept>
#include <assert.h>
#include <iostream>
#include <fstream>
#include <algorithm>  

using namespace std;

namespace Terran {

MethodsPeriodicGaussian::MethodsPeriodicGaussian(const vector<Param> &params, 
    double period) : Methods(params), period_(period) {

}

// Newton-Raphson using mean of the components to find the maxima (and minima)
/*
vector<double> MethodsGaussian::findMaxima() const {
    vector<double> maxima;

    double min =  numeric_limits<double>::max();
    double max = -numeric_limits<double>::max();
    double max_sig = -numeric_limits<double>::max();
    for(int i=0; i < params_.size(); i++) {
        min = (params_[i].u < min) ? params_[i].u : min;
        max = (params_[i].u > max) ? params_[i].u : max;
        max_sig = (params_[i].s > max) ? params_[i].s : max;
    }

    while(gaussianMixture(params_,min) > 1e-5) 
        min -= max_sig;
    while(gaussianMixture(params_,max) > 1e-5) 
        max += max_sig;

    for(int k=0; k<params_.size(); k++) {
        double x_old = params_[k].u;
        bool found = false;
        int iteration = 0;
        while(!found) {
            iteration++;
            if(iteration >= 1e3) {
                throw(std::runtime_error("findMaxima: maximum iteration count reached!"));
            }
            double x_new = x_old - gaussianMixtureDx(params_, x_old)/gaussianMixtureDx2(params_, x_old);
            if(fabs(x_new-x_old) < 1e-8)
                found = true;
            if(x_new < min || x_new > max)
                break;
            x_old = x_new;
        }

        // see if found point is a maxima
        if(gaussianMixtureDx(params_, x_old-delta) > 0 && gaussianMixtureDx(params_, x_old+delta) < 0) {
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
*/

/*
vector<double> MethodsPeriodicGaussian::findMaxima() const {
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
                throw(std::runtime_error("findMaxima: maximum iteration count reached!"));
            }
            double xn_new = xn_old + delta*periodicGaussianMixtureDx(params_, xn_old, period_);
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
*/

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
                throw(std::runtime_error("findMaxima: maximum iteration count reached!"));
            }
            double x_new = x_old - periodicGaussianMixtureDx(params_, x_old, period_)/periodicGaussianMixtureDx2(params_, x_old, period_);
            
            x_new = normalize(x_new);
            
            if(fabsp(x_new, x_old, period_) < 1e-8)
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
        //} while (fabs(my) > 1e-9);
        } while (fabsp(ak,bk,period_) > 1e-8);
        minima.push_back(mk);
    }
    return minima;
}

}
