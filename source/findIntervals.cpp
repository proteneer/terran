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

// Maximas can be found by using the mean as initial guesses. In the case of overlapping components
// They should converge to the same point. Unlike minimas, they don't have pathological zero derivatives
vector<double> findPeriodicMaxima(const vector<Param> &params, double period, int images) {
    vector<double> maximas;
    for(int k=0; k<params.size(); k++) {
        double xn_old = params[k].u;
        bool found = false;
        const double delta = 1e-4;
        int iteration = 0;
        while(!found) {
            iteration++;
            // Saddle points suck. If the simple algorithm doesn't converge
            // it's probably due to a saddle point. 
            if(iteration >= 1e7) {
                throw(std::runtime_error("findPeriodicMaxima: maximum iteration count reached!"));
            }

            double xn_new = xn_old + delta*periodicGaussianMixtureDx(params, xn_old, period, images);
            if(fabs(xn_new - xn_old) < 1e-8)
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

static double normalize(double x, double left=-PI/2, double right=PI/2) {
    double period = right-left;
    while(x > right) {
        x -= period;
    }
    while(x < left) {
        x += period;
    }
    return x;
}

// Algorithm 1.

// 1. Find all the maximas
// 2. Sort the maximas according to their x value
// 3. Find minimas by taking adjacent pair of maximas and doing either
//    bisection or gradient descent to converge at a minima
//
// Note: Be careful of corner cases dealing with periodic boundaries. 
vector<double> findPeriodicMinima(const vector<Param> &params, double period, int images) {

    // Get a list of all the maxima and sort them
    vector<double> maxima = findPeriodicMaxima(params, period, images);

    assert(maxima.size() > 1);

    sort(maxima.begin(), maxima.end());

    for(int i=0; i < maxima.size()-1; i++) {

    }

    // Take care of the corner case with the boundary maxima
    double l1;
    {
        // perturbation so they fall down the right direction
        double pt = 0.01;
        double xn_old = maxima[0];
        xn_old += pt;
        xn_old = normalize(xn_old);
        double delta = 1e-4;
        int iteration = 0;
        bool found = false;
        // gradient descent from the xn_old point
        while(!found) {
            iteration++;
            // Saddle points suck. If the simple algorithm doesn't converge
            // it's probably due to a saddle point. 
            if(iteration >= 1e10) {
                //throw(std::runtime_error("findPeriodicMinima Left: maximum iteration count reached!"));
            }
            double xn_new = xn_old - delta*periodicGaussianMixtureDx(params, xn_old, period, images);
            if(fabs(xn_new - xn_old) < 1e-8)
                found = true;
            xn_old = xn_new;
        }
        l1 = xn_old;
    }
    double r1;
    {
        double pt = 0.01;
        double xn_old = maxima[maxima.size()-1];
        xn_old -= pt;
        xn_old = normalize(xn_old);
        double delta = 1e-4;
        int iteration = 0;
        bool found = false;
        // gradient descent from the xn_old point
        while(!found) {
            iteration++;
            // Saddle points suck. If the simple algorithm doesn't converge
            // it's probably due to a saddle point. 
            if(iteration >= 1e10) {
                //throw(std::runtime_error("findPeriodicMaxima Right: maximum iteration count reached!"));
            }
            double xn_new = xn_old - delta*periodicGaussianMixtureDx(params, xn_old, period, images);
            if(fabs(xn_new - xn_old) < 1e-8)
                found = true;
            xn_old = xn_new;
        }
        r1 = xn_old;
    }

    double min = normalize(((l1+period)+r1)/2);

    cout << l1 << " " << r1 << " " << min << endl;

    // Eg:
    // Given maximas, a,b,c,d,e
    // There must be minimas between:
    // (e,a) (a,b) (c,d) (d,e)
    return vector<double>(0);
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
