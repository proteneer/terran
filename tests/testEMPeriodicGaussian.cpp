// tests expectation maximization of periodic gaussians for correctness

#include <math.h>
#include <vector>
#include <iostream>
#include <fstream>
#include <algorithm>

#include <MethodsPeriodicGaussian.h>
#include <EMPeriodicGaussian.h>
#include <MathFunctions.h>
#include <sstream>

#include "util.h"

using namespace std;
using namespace Terran;

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
    em.run();
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
    for(int i=0; i<5000*initParams[0].p; i++) {
        data.push_back(periodicGaussianSample(initParams[0].u,initParams[0].s,period));
    }
    initParams[1].p = 0.35;
    initParams[1].u = 1.9;
    initParams[1].s = 0.5;
    for(int i=0; i<5000*initParams[1].p; i++) {
        data.push_back(periodicGaussianSample(initParams[1].u,initParams[1].s,period));
    }

    random_shuffle(data.begin(), data.end());
    data.resize(2500);

    // the rough range of the maxima and minima of this distribution using the histogram
    vector<double> truthMaxima, truthMaximaErrors;
    truthMaxima.push_back(-0.245); // += 0.12
    truthMaximaErrors.push_back(0.12);
    truthMaxima.push_back( 1.895); // += 0.21
    truthMaximaErrors.push_back(0.21);
    vector<double> truthMinima, truthMinimaErrors;
    truthMinima.push_back(-2.38); // += 0.37
    truthMinimaErrors.push_back(0.37);
    truthMinima.push_back( 0.99); // += 0.31
    truthMinimaErrors.push_back(0.31);

    vector<double> maxima, minima;
    EMPeriodicGaussian em(data, period);


    try {
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
        em.setParameters(params);
        em.setTolerance(0.01);
        em.run();
        vector<Param> optimizedParams = em.getParams();
        Util::matchParameters(initParams, optimizedParams, 0.05);
        MethodsPeriodicGaussian mpg(optimizedParams, period);
        maxima = mpg.findMaxima();
        minima = mpg.findMinima();
        Util::matchPoints(truthMaxima, maxima, truthMaximaErrors);
        Util::matchPoints(truthMaxima, maxima, truthMaximaErrors);
    } catch(const exception &e) {
        cout << "bimodal run() test failed!" << endl;
        throw e;
    }


    try{
        vector<Param> params;
        const int numSamples = 8;
        for(int i=0; i < numSamples; i++) {
            Param p;
            p.p = 1/(float)numSamples;
            p.u = -PI + ((float)i/numSamples)*(period);
            p.s = 0.3;
            params.push_back(p);
        }
        em.setParameters(params);
        em.simpleRun(25);
        vector<Param> optimizedParams = em.getParams();
        MethodsPeriodicGaussian mpg(optimizedParams, period);
        maxima = mpg.findMaxima();
        minima = mpg.findMinima();
        Util::matchPoints(truthMaxima, maxima, truthMaximaErrors);
        Util::matchPoints(truthMaxima, maxima, truthMaximaErrors);
    } catch(const exception &e) {
        cout << "bimodal simpleRun() test failed!" << endl;
        throw e;
    }
}

int main() {
    try {
        srand(1);
        testUnimodalPeriodicGaussian();
        srand(1);
        testBimodalPeriodicGaussian();
    } catch( const std::exception &e ) {
        cout << e.what() << endl;
    }
    cout << "done" << endl;
}
