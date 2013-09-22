#include "MethodsGaussian.h"
#include <stdexcept>
#include <assert.h>
#include <iostream>
#include <fstream>
#include <algorithm>  

using namespace std;

namespace Terran {

// useful utility for plotting the gaussin
static void plotGaussian(const vector<Param> &params) {
    ofstream mixture("mixture.dat");
    for(double xn = -10; xn < 10; xn += 0.01) {
        mixture << xn << " " << gaussianMixture(params, xn) << endl;
    }
    ofstream mixtureDx("mixtureDx.dat");
    for(double xn = -10; xn < 10; xn += 0.01) {
        mixtureDx << xn << " " << gaussianMixtureDx(params, xn) << endl;
    }
}

MethodsGaussian::MethodsGaussian(const vector<Param> &params) : Methods(params) {

}

// Newton-Raphson using mean of the components to find the maxima (and minima)
vector<double> MethodsGaussian::findMaxima() const {
    vector<double> maximas;

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
        const double delta = 1e-4;
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
            for(int i=0; i<maximas.size(); i++)
                if(fabs(x_old-maximas[i]) < 1e-3)
                    skip = true;
            if(!skip) 
                maximas.push_back(x_old);
        }
    }
    return maximas;
}

// TODO: This does redundant work, can adapt findMaxima() to do this.
vector<double> MethodsGaussian::findMinima() const {
    vector<double> minima;
    vector<double> maxima = findMaxima();
    assert(maxima.size() > 0);
    sort(maxima.begin(), maxima.end());
    const double perturb = 5e-2;

    if(maxima.size() > 1) {
        for(int k=0; k<maxima.size()-1; k++) {
            double ak = maxima[k];
            double bk = maxima[k+1];
            ak += perturb;
            bk -= perturb;
            double mk = 0;
            double my = 0;
            int iteration = 0;
            // bisect the interval
            do {
                if(iteration > 1e5) {
                    throw(std::runtime_error("MethodsGaussian::findMinima() - maximized number of iterations reached."));
                }
                else
                    iteration++;
                
                mk = (ak+bk)/2.0;
                my = gaussianMixtureDx(params_, mk);
                if(my < 0) 
                    ak = mk;
                else
                    bk = mk;
            } while (fabs(ak-bk) > 1e-8);
            minima.push_back(mk);
        }
    }
    return minima;
}

} // namespace Terran
