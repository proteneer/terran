// tests expectation maximization of gaussians for correctness

#include <math.h>
#include <vector>
#include <iostream>
#include <fstream>
#include <sstream>
#include <algorithm>

#include <EMGaussian.h>
#include <EMPeriodicGaussian.h>
#include <MathFunctions.h>
#include <MethodsGaussian.h>

#include "util.h"

using namespace std;
using namespace Terran;

// returns a gaussian sample using the Box-Muller transform
// where, x1,x2 ~ U[0,1)
double2 boxMullerSample(double x1, double x2, double u, double s) {
    double2 z;
    z.x = sqrt(-2*log(x1))*cos(2*PI*x2);
    z.y = sqrt(-2*log(x1))*sin(2*PI*x2);
    z.x *= s;
    z.x += u;
    z.y *= s;
    z.y += u;
    return z;
}

void testUnimodalGaussian() {
    vector<Param> initParams(1);
    initParams[0].p = 1;
    initParams[0].u = 12.3;
    initParams[0].s = 4.7;
    vector<double> data;
    for(int i=0; i<10000; i++) {
        double x1 = (double)rand()/(double)RAND_MAX;
        double x2 = (double)rand()/(double)RAND_MAX;
        double2 z = boxMullerSample(x1, x2, initParams[0].u, initParams[0].s);
        data.push_back(z.x);
        data.push_back(z.y);
    }
    vector<Param> params;
    Param p;
    p.p = 1;
    p.u = 6.2;
    p.s = 8.1;
    params.push_back(p);
    EMGaussian em(data, params);
    em.run();
    vector<Param> optimizedParams = em.getParams();
    Util::matchParameters(initParams, optimizedParams, 0.1);
}

void testBimodalGaussian() {
    int numSamples = 20000;
    vector<double> data;
    vector<Param> initParams(2);
    initParams[0].p = 0.4;
    initParams[0].u = -3.4;
    initParams[0].s = 1.2;
    for(int i=0; i<numSamples*initParams[0].p; i++) {
        double x1 = (double)rand()/(double)RAND_MAX;
        double x2 = (double)rand()/(double)RAND_MAX;
        double2 z = boxMullerSample(x1, x2, initParams[0].u, initParams[0].s);
        data.push_back(z.x);
        data.push_back(z.y);
    }

    initParams[1].p = 0.6;
    initParams[1].u = 7.4;
    initParams[1].s = 6.2;
    for(int i=0; i<numSamples*initParams[1].p; i++) {
        double x1 = (double)rand()/(double)RAND_MAX;
        double x2 = (double)rand()/(double)RAND_MAX;
        double2 z = boxMullerSample(x1, x2, initParams[1].u, initParams[1].s);
        data.push_back(z.x);
        data.push_back(z.y);
    }
    vector<Param> params;
    {
        Param p;
        p.p = 0.5;
        p.u = 6.2;
        p.s = 8.1;
        params.push_back(p);
    }
    {
        Param p;
        p.p = 0.5;
        p.u = 0;
        p.s = 2.1;
        params.push_back(p);
    }
    EMGaussian em(data, params);
    em.run();
    vector<Param> optimizedParams = em.getParams();
    Util::matchParameters(initParams, optimizedParams, 0.1);
}

void testSimpleRunBimodal() {
    int numSamples = 20000;
    vector<double> data;
    vector<Param> initParams(2);
    initParams[0].p = 0.4;
    initParams[0].u = -3.4;
    initParams[0].s = 1.2;
    for(int i=0; i<numSamples*initParams[0].p; i++) {
        double x1 = (double)rand()/(double)RAND_MAX;
        double x2 = (double)rand()/(double)RAND_MAX;
        double2 z = boxMullerSample(x1, x2, initParams[0].u, initParams[0].s);
        data.push_back(z.x);
        data.push_back(z.y);
    }

    initParams[1].p = 0.6;
    initParams[1].u = 7.4;
    initParams[1].s = 6.2;
    for(int i=0; i<numSamples*initParams[1].p; i++) {
        double x1 = (double)rand()/(double)RAND_MAX;
        double x2 = (double)rand()/(double)RAND_MAX;
        double2 z = boxMullerSample(x1, x2, initParams[1].u, initParams[1].s);
        data.push_back(z.x);
        data.push_back(z.y);
    }
    vector<Param> params;
    {
        Param p;
        p.p = 0.5;
        p.u = 6.2;
        p.s = 8.1;
        params.push_back(p);
    }
    {
        Param p;
        p.p = 0.5;
        p.u = 0;
        p.s = 2.1;
        params.push_back(p);
    }
    vector<double> subset(data);
    random_shuffle(subset.begin(), subset.end());
    subset.resize(2500);
    EMGaussian em(subset, params);

    bool rc = em.simpleRun(25);

    vector<Param> optimizedParams = em.getParams();
    MethodsGaussian method(optimizedParams);
    vector<double> maxima = method.findMaxima();
    vector<double> minima = method.findMinima();

    vector<double> trueMaxima;
    trueMaxima.push_back(-3.25);
    trueMaxima.push_back( 7.4);
    vector<double> errorsMaxima;
    errorsMaxima.push_back(0.5);
    errorsMaxima.push_back(1.5);
    Util::matchPoints(maxima, trueMaxima, errorsMaxima);
    vector<double> trueMinima;
    trueMinima.push_back(0.5);
    Util::matchPoints(minima, trueMinima, 0.3);

    if(rc != 1) {
        throw std::runtime_error("EM - Failed to converge");
    }    
}

int main() {
    try {
		testSimpleRunBimodal();
		srand(1);
        testUnimodalGaussian();
        srand(1);
        testBimodalGaussian();
    } catch( const std::exception &e ) {
        cout << e.what() << endl;
    }
    cout << "done" << endl;
}
