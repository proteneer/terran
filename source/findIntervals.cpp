#include "findIntervals.h"

#include <vector>

using namespace std;

void partitionGaussian(const vector<Param> &params) {


}


#include <iostream>


void partitionPeriodicGaussian(const vector<Param> &params, double period, int images) {

/*
    for(double xn = -period/2; xn < period/2; xn += 0.01) {
        cout << xn << " " << periodicGaussianMixture(params, xn, period, images) << endl;
    }
*/
    for(double xn = -period/2; xn < period/2; xn += 0.01) {
        cout << xn << " " << periodicGaussianMixtureDx(params, xn, period, images) << endl;
    }

}
