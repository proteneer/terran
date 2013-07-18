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

// Golden section method of finding minima.
// Original implementation: Numerical Recipes
// Modified for periodic boundaries
vector<double> findPeriodicMinimaGS(const vector<Param> &params, double period, int images) {
    vector<double> minima;
    vector<double> maxima = findPeriodicMaxima(params, period, images);
    assert(maxima.size() > 0);
    sort(maxima.begin(), maxima.end());

    const double tol = 1e-6;
    const double R = 0.61803399;
    const double C = (1.0-R);

    for(int i=0; i < maxima.size(); i++) {
        // x0 will always be ax
        // x3 will always be cx
        // x1 and x2 represent bx and the guess
        double x0,x1,x2,x3;
        double ax = maxima[i];
        double cx = maxima[(i+1)%maxima.size()];
        x0 = ax;
        x3 = cx;

        if(i == maxima.size()-1) {
            ax += period;
        } 
        double bx = (ax+cx)/2;
        bx = normalize(bx); 
        cout << "a b c: " << ax << " " << bx << " " << cx << endl;
        if(fabsp(cx,bx,period) > fabsp(bx,ax,period)) {
            x1 = bx;
            x2 = bx + C*(cx-bx);
        } else {
            x2 = bx;
            x1 = bx - C*(bx-ax);
        }
        // move into periodic box
        x1 = normalize(x1);
        x2 = normalize(x2);
        double f1 = periodicGaussianMixture(params, x1, period, images);
        double f2 = periodicGaussianMixture(params, x2, period, images);
        cout << "Initial: " << x0 << " " << x1 << " " << x2 << " " << x3 << endl;

        // do loop is critical here, else the first minima unit test will fail!
        do {
            cout << std::setprecision(8) << x0 << " " << x1 << " " << x2 << " " << x3 << endl;
            if(f2 < f1) {
                x0 = x1;
                x1 = x2;
                x2 = R*x1+C*x3;
                normalize(x2);
                f1 = f2;
                f2 = periodicGaussianMixture(params, x2, period, images);
            } else {
                x3 = x2;
                x2 = x1;
                x1 = R*x2+C*x0;
                normalize(x1);
                f2 = f1;
                f1 = periodicGaussianMixture(params, x1, period, images);
            }
        } while(fabsp(x3,x0,period) > tol);
        minima.push_back(x3);
    }
    return minima;
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
    vector<double> minima;
    vector<double> maxima = findPeriodicMaxima(params, period, images);
    assert(maxima.size() > 0);
    sort(maxima.begin(), maxima.end());
    for(int i=0; i < maxima.size(); i++) {
        double l1;
        const unsigned long maxIterations = 1e11;
        {
            // add small perturbation
            double pt = 0.01;
            double xn_old = maxima[i];
            xn_old = normalize(xn_old+pt);
            double delta = 1e-4;
            unsigned long iteration = 0;
            bool found = false;
            // gradient descent from the xn_old point
            while(!found) {
                iteration++;
                if(iteration >= maxIterations) {
                    //throw(std::runtime_error("findPeriodicMinima Left: maximum iteration count reached!"));
                }
                double xn_new = xn_old - delta*periodicGaussianMixtureDx(params, xn_old, period, images);
                xn_new = normalize(xn_new);
                if(fabsp(xn_new,xn_old,period) < 1e-8)
                    found = true;
                xn_old = xn_new;
            }
            cout << "left converged in " << iteration << endl;
            l1 = xn_old;
        }
        double r1;
        {
            double pt = 0.01;
            double xn_old = maxima[(i+1)%maxima.size()];
            xn_old = normalize(xn_old-pt);
            double delta = 1e-4;
            unsigned long iteration = 0;
            bool found = false;
            // gradient descent from the xn_old point
            while(!found) {
                iteration++;
                if(iteration >= maxIterations) {
                    //throw(std::runtime_error("findPeriodicMaxima Right: maximum iteration count reached!"));
                }
                double xn_new = xn_old - delta*periodicGaussianMixtureDx(params, xn_old, period, images);
                xn_new = normalize(xn_new);
                if(fabsp(xn_new,xn_old,period) < 1e-8)
                    found = true;
                xn_old = xn_new;
            }
            cout << "right converged in " << iteration << endl;
            r1 = xn_old;
        }
        double min = normalize(l1 + fabsp(l1,r1,period)/2);
        //cout << l1 << " " << r1 << " " << min << endl;
        minima.push_back(min);
    }
    // Eg:
    // Given maximas, a,b,c,d,e
    // There must be minimas between:
    // (e,a) (a,b) (c,d) (d,e)
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
