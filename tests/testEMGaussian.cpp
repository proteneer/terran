// tests expectation maximization of gaussians for correctness

#include <math.h>
#include <vector>
#include <iostream>

#include <EMGaussian.h>
#include <MathFunctions.h>

#include "util.h"

using namespace std;

struct double2 {
    double x;
    double y;
};

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
    em.run(10000, 0.1);
    vector<Param> optimizedParams = em.getParams();
    Util::matchParameters(initParams, optimizedParams, 0.02);
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
    em.run(10000, 0.1);
    vector<Param> optimizedParams = em.getParams();
    Util::matchParameters(initParams, optimizedParams, 0.05);
}

void testTrimodalGaussian() {

}

int main() {
    try {
        testUnimodalGaussian();
        testBimodalGaussian();
    } catch( const std::exception &e ) {
        cout << e.what() << endl;
    }
    cout << "done" << endl;
}
