// tests expectation maximization of gaussians for correctness

#include <math.h>
#include <vector>
#include <iostream>

#include <EMGaussian.h>
#include <MathFunctions.h>

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
    vector<double> data;
    double u1 = 12.3;
    double s1 = 4.7;
    for(int i=0; i<10000; i++) {
        double x1 = (double)rand()/(double)RAND_MAX;
        double x2 = (double)rand()/(double)RAND_MAX;
        double2 z = boxMullerSample(x1, x2, u1, s1);
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
    if(fabs(optimizedParams[0].u - u1) > 0.2) {
        throw(std::runtime_error("testUnimodalGaussian failed, u mismatch"));
    }
    if(fabs(optimizedParams[0].s - s1) > 0.2) {
        throw(std::runtime_error("testUnimodalGaussian failed, s mismatch"));
    }
}

void testBimodalGaussian() {
    vector<double> data;
    double u1 = -3.4;
    double s1 = 1.2;
    for(int i=0; i<10000; i++) {
        double x1 = (double)rand()/(double)RAND_MAX;
        double x2 = (double)rand()/(double)RAND_MAX;
        double2 z = boxMullerSample(x1, x2, u1, s1);
        data.push_back(z.x);
        data.push_back(z.y);
    }
    double u2 = 7.4;
    double s2 = 6.2;
    for(int i=0; i<20000; i++) {
        double x1 = (double)rand()/(double)RAND_MAX;
        double x2 = (double)rand()/(double)RAND_MAX;
        double2 z = boxMullerSample(x1, x2, u2, s2);
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
    if(fabs(optimizedParams[0].u - u1) > 0.2 ) {
        throw(std::runtime_error("testUnimodalGaussian failed, u mismatch"));
    }
    if(fabs(optimizedParams[0].s - s1) > 0.2) {
        throw(std::runtime_error("testUnimodalGaussian failed, s mismatch"));
    }
    if(fabs(optimizedParams[1].u - u1) > 0.2 ) {
        throw(std::runtime_error("testUnimodalGaussian failed, u mismatch"));
    }
    if(fabs(optimizedParams[1].s - s1) > 0.2) {
        throw(std::runtime_error("testUnimodalGaussian failed, s mismatch"));
    }
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
