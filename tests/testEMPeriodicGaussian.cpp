// tests expectation maximization of periodic gaussians for correctness

#include <math.h>
#include <vector>
#include <iostream>

#include <EMPeriodicGaussian.h>
#include <MathFunctions.h>

#include <algorithm>

using namespace std;

// draws a random sample from the periodic gaussian with mean u,
// std. deviation s, and period period
// assume x is in -period/2, period/2
double periodicGaussianSample(double u, double s, double period) {
    // sample between [0,1]
    while(true) {
        double x1 = (double) rand() / (double) RAND_MAX;
        // move to [-0.5, 0.5]
        x1 -= 0.5;
        // move to [-0.5*2PI, -0.5*2PI]
        x1 *= period;
        // determine the probability of drawing x1
        double probability = periodicGaussian(u, s, x1, 100, period);
        // make a random sample between 0 and 1 
        double draw = (double) rand() / (double) RAND_MAX;
        if(draw < probability) {
            return x1;
        }
    }
}

void testUnimodalPeriodicGaussian() {
    vector<double> data;
    double period = 2*PI;
    double u1 = 1.2345;
    double s1 = 1.1;
    for(int i=0; i < 5000; i++) {
        data.push_back(periodicGaussianSample(u1,s1, period));
    }

    sort(data.begin(), data.end());

    Param p;
    p.p = 1;
    p.u = 1.7;
    p.s = 2.1;
    vector<Param> params;
    params.push_back(p);
    EMPeriodicGaussian em(data, params);
    em.EStep();
    em.MStep();
}

int main() {
    try {
        testUnimodalPeriodicGaussian();
    } catch( const std::exception &e ) {
        cout << e.what() << endl;
    }
    cout << "done" << endl;
}
