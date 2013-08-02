// tests expectation maximization of periodic gaussians for correctness

#include <math.h>
#include <vector>
#include <iostream>
#include <fstream>

#include <MethodsPeriodicGaussian.h>
#include <EMPeriodicGaussian.h>
#include <MathFunctions.h>

#include "util.h"

using namespace std;
using namespace Terran;

void testSampleDomain() {
    double period = 2*PI;
    vector<double> data;
    data.push_back(0);


}

void testUnimodalPeriodicGaussian() {
    vector<double> data;
    vector<Param> initParams(1);
    double period = 2*PI;
    initParams[0].p = 1;
    initParams[0].u = 1.2345;
    initParams[0].s = 1.1;
    for(int i=0; i < 5000; i++) {
        data.push_back(periodicGaussianSample(initParams[0].u,initParams[0].s,period));
    }
    Param p;
    p.p = 1;
    p.u = 1.7;
    p.s = 2.1;
    vector<Param> params;
    params.push_back(p);
    EMPeriodicGaussian em(data, params, period);
    em.run(100, 0.1);
    vector<Param> optimizedParams = em.getParams();
    Util::matchParameters(initParams, optimizedParams, 0.05);
}

void testBimodalPeriodicGaussian() {
    vector<double> data;
    double period = 2*PI;
    vector<Param> initParams(2);
    initParams[0].p = 0.65;
    initParams[0].u = -0.3;
    initParams[0].s = 0.5;
    for(int i=0; i<10000*initParams[0].p; i++) {
        data.push_back(periodicGaussianSample(initParams[0].u,initParams[0].s,period));
    }
    initParams[1].p = 0.35;
    initParams[1].u = 1.9;
    initParams[1].s = 0.5;
    for(int i=0; i<10000*initParams[1].p; i++) {
        data.push_back(periodicGaussianSample(initParams[1].u,initParams[1].s,period));
    }
    
    vector<Param> params;
    {
        Param p;
        p.p = 0.5;
        p.u = -1.0;
        p.s = 0.5;
        params.push_back(p);
    }
    {
        Param p;
        p.p = 0.5;
        p.u = 2.1;
        p.s = 1.4;
        params.push_back(p);
    }
    EMPeriodicGaussian em(data, params, period);
    em.run(100, 0.1);
    vector<Param> optimizedParams = em.getParams();
    Util::plotPeriodicGaussian(initParams, period, "bimodal");
    Util::matchParameters(initParams, optimizedParams, 0.05);
}


void testOverfitPeriodicGaussian() {
    vector<double> data;
    double period = 2*PI;
    vector<Param> initParams(2);
    initParams[0].p = 0.65;
    initParams[0].u = -0.3;
    initParams[0].s = 0.5;
    for(int i=0; i<5000*initParams[0].p; i++) {
        data.push_back(periodicGaussianSample(initParams[0].u,initParams[0].s,period));
    }
    initParams[1].p = 0.35;
    initParams[1].u = 1.9;
    initParams[1].s = 0.5;
    for(int i=0; i<5000*initParams[1].p; i++) {
        data.push_back(periodicGaussianSample(initParams[1].u,initParams[1].s,period));
    }
    
    /*
    ofstream d("data.txt");
    for(int i=0; i<data.size(); i++) {
         d << data[i] << endl;
    }
    */

    vector<Param> params;
    const int numSamples = 8;
    for(int i=0; i < numSamples; i++) {
        Param p;
        p.p = 1/(float)numSamples;
        p.u = -PI + ((float)i/numSamples)*(period);
        p.s = 0.3;
        params.push_back(p);
    }

    EMPeriodicGaussian em(data, params, period);
    em.run(1000,0.001);
    vector<Param> optimizedParams = em.getParams();

    MethodsPeriodicGaussian mpg(optimizedParams, period);
    vector<double> maxima = mpg.findMaxima();
    vector<double> minima = mpg.findMinima();

    vector<double> truthMaxima;
    truthMaxima.push_back( 1.898);
    truthMaxima.push_back(-0.301);
    vector<double> truthMinima;
    truthMinima.push_back( 0.888);
    truthMinima.push_back(-2.382);

    double tol = period/100.0;

    Util::matchPoints(truthMaxima, maxima, tol);
    Util::matchPoints(truthMinima, minima, tol);
}

void testAdaptiveRun() {
    vector<double> data;
    double period = 2*PI;
    vector<Param> initParams(2);
    initParams[0].p = 0.65;
    initParams[0].u = -0.3;
    initParams[0].s = 0.5;
    for(int i=0; i<5000*initParams[0].p; i++) {
        data.push_back(periodicGaussianSample(initParams[0].u,initParams[0].s,period));
    }
    initParams[1].p = 0.35;
    initParams[1].u = 1.9;
    initParams[1].s = 0.5;
    for(int i=0; i<5000*initParams[1].p; i++) {
        data.push_back(periodicGaussianSample(initParams[1].u,initParams[1].s,period));
    }
    
    ofstream d("data.txt");
    for(int i=0; i<data.size(); i++) {
         d << data[i] << endl;
    }

    vector<Param> params;
    const int numSamples = 8;
    for(int i=0; i < numSamples; i++) {
        Param p;
        p.p = 1/(float)numSamples;
        p.u = -PI + ((float)i/numSamples)*(period);
        p.s = 0.3;
        params.push_back(p);
    }

    EMPeriodicGaussian em(data, params, period);
    em.adaptiveRun(1000,0.01, 0.08);
    vector<Param> optimizedParams = em.getParams();

    MethodsPeriodicGaussian mpg(optimizedParams, period);
    vector<double> maxima = mpg.findMaxima();
    vector<double> minima = mpg.findMinima();

    Util::plotPeriodicGaussian(optimizedParams, period, "adaptive");

    vector<double> truthMaxima;
    truthMaxima.push_back( 1.898);
    truthMaxima.push_back(-0.301);
    vector<double> truthMinima;
    truthMinima.push_back( 0.888);
    truthMinima.push_back(-2.382);

/*
    cout << "maxima" << endl;
    for(int i=0; i < maxima.size(); i++) {
        cout << maxima[i] << endl;
    }

    cout << "minima" << endl;
    for(int i=0; i < minima.size(); i++) {
        cout << minima[i] << endl;
    }
*/

    double tol1 = period/100.0f;
    double tol2 = period/50.0f;

    Util::matchPoints(truthMaxima, maxima, tol1);
    Util::matchPoints(truthMinima, minima, tol2);
}

int main() {
    try {
        testSampleDomain();
        srand(1);
        testUnimodalPeriodicGaussian();
        srand(1);
        testBimodalPeriodicGaussian();
        srand(1);
        //testOverfitPeriodicGaussian();
        srand(1);
        testAdaptiveRun();
    } catch( const std::exception &e ) {
        cout << e.what() << endl;
    }
    cout << "done" << endl;
}
