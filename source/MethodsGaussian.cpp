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
    ofstream mixtureDx2("mixtureDx2.dat");
    for(double xn = -10; xn < 10; xn += 0.01) {
        mixtureDx2 << xn << " " << gaussianMixtureDx2(params, xn) << endl;
    }
}

MethodsGaussian::MethodsGaussian(const vector<Param> &params) : Methods(params) {

}

// Newton-Raphson doesn't work - no guarantee of global convergence
vector<double> MethodsGaussian::findMaxima() const {
 vector<double> maximas;
    for(int k=0; k<params_.size(); k++) {
        double xn_old = params_[k].u;
        bool found = false;
        const double delta = 1e-2;
        int iteration = 0;
        while(!found) {
            iteration++;
            // Saddle points suck. If the simple algorithm doesn't converge
            // it's probably due to a saddle point. 
            if(iteration >= 1e5) {
                throw(std::runtime_error("findMaxima: maximum iteration count reached!"));
            }
            double xn_new = xn_old + delta*gaussianMixtureDx(params_, xn_old);
            if(fabs(gaussianMixtureDx(params_, xn_new)) < 1e-4)
                found = true;
            xn_old = xn_new;
        }
        bool skip = false;
        // if two gaussians are too close to each other, discard it
        for(int i=0; i<maximas.size(); i++)
            if(fabs(xn_old-maximas[i]) < 1e-3)
                skip = true;
        if(!skip) 
            maximas.push_back(xn_old);
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
