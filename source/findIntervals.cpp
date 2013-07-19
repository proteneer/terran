#include "findIntervals.h"
#include <iomanip>

#include <vector>

#include <iostream>
#include <fstream>
#include <stdexcept>

#include <assert.h>
using namespace std;

void partitionGaussian(const vector<Param> &params) {

}

// Maxima can be found by using the mean as initial guesses. In the case of overlapping components
// They should converge to the same point. Unlike minimas, they don't have pathological zero derivatives
vector<double> findPeriodicMaxima(const vector<Param> &params, double period, int images) {
    vector<double> maximas;
    for(int k=0; k<params.size(); k++) {
        double xn_old = params[k].u;
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
            double xn_new = xn_old + delta*periodicGaussianMixtureDx(params, xn_old, period, images);
            xn_new = normalize(xn_new);
            if(fabsp(xn_new, xn_old, period) < 1e-8)
                found = true;
            xn_old = xn_new;
        }
        bool skip = false;
        // if two gaussians are too close to each other, discard it
        for(int i=0; i<maximas.size(); i++)
            if(fabsp(xn_old, maximas[i], period) < 1e-3)
                skip = true;
        if(!skip) 
            maximas.push_back(xn_old);
    }
    return maximas;
}

// Find the minima via the bisection of the periodicGaussianMixture derivative
// So far, this is the best method, superior to that of Golden Section Search and Gradient Descent
vector<double> findPeriodicMinimaBS(const vector<Param> &params, double period, int images) {
    vector<double> minima;
    vector<double> maxima = findPeriodicMaxima(params, period, images);
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
                mk = normalize((ak+period+bk)/2.0);
            else
                mk = normalize((ak+bk)/2.0);
            my = periodicGaussianMixtureDx(params, mk, period, images);
            if(my < 0) 
                ak = mk;
            else
                bk = mk;
        //} while (fabs(my) > 1e-9);
        } while (fabsp(ak,bk,period) > 1e-8);
        minima.push_back(mk);
    }
    return minima;
}

void partitionPeriodicGaussian(const vector<Param> &params, double period, int images) {
    ofstream mixture("mixture.dat");
    for(double xn = -period/2; xn < period/2; xn += 0.01) {
        mixture << xn << " " << periodicGaussianMixture(params, xn, period, images) << endl;
    }
    ofstream mixtureDx("mixtureDx.dat");
    for(double xn = -period/2; xn < period/2; xn += 0.01) {
        mixtureDx << xn << " " << periodicGaussianMixtureDx(params, xn, period, images) << endl;
    }
    /*
    ofstream mixture("mixture.dat");
    double left = 1.4732;
    double right = 1.4734;
    for(double xn = left; xn < right; xn += 0.00001) {
        mixture << std::setprecision(10) << xn << " " << periodicGaussianMixture(params, xn, period, images) << endl;
    }
    ofstream mixtureDx("mixtureDx.dat");
    for(double xn = left; xn < right; xn += 0.00001) {
        mixtureDx << std::setprecision(10) << xn << " " << periodicGaussianMixtureDx(params, xn, period, images) << endl;
    }
    */
}
