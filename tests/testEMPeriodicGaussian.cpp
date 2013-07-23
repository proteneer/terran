// tests expectation maximization of periodic gaussians for correctness

#include <math.h>
#include <vector>
#include <iostream>

#include <EMPeriodicGaussian.h>
#include <MathFunctions.h>

#include <algorithm>

#include "util.h"

using namespace std;

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
    EMPeriodicGaussian em(data, params);
    em.run(10000, 0.1);
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
    
    EMPeriodicGaussian em(data, params, period, 5);
    em.run(10000, 0.1);
    vector<Param> optimizedParams = em.getParams();

    Util::matchParameters(initParams, optimizedParams, 0.05);
}

int main() {
    try {
        testUnimodalPeriodicGaussian();
        testBimodalPeriodicGaussian();
    } catch( const std::exception &e ) {
        cout << e.what() << endl;
    }
    cout << "done" << endl;
}
