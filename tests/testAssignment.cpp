#include <math.h>
#include <vector>
#include <iostream>

#include <EMGaussian.h>
#include <EMPeriodicGaussian.h>
#include <PartitionPeriodicGaussian.h>
#include <MathFunctions.h>

#include "util.h"

using namespace std;

// A three-cluster dataset looking something like this:

// PI
// |       *
// |      ***
// |     **0**
// |      ***     *
// |       *     ***
// |            **2**
// |       *     ***
// |      ***     *
// |     **1**
// |      ***
// |       *
//-PI-------------------PI


void testEasyCase2D() {

    vector<vector<double> > dataset;
    // setup 0th cluster;
    {
        double u1 = -PI/2;
        double u2 =  PI/2;
        double s = 0.3;
        double period = 2*PI;
        for(int i=0; i < 4000; i++) {
            vector<double> point(2);
            point[0] = periodicGaussianSample(u1, s, period);
            point[1] = periodicGaussianSample(u2, s, period);
            dataset.push_back(point);
        }
    }

    // setup 1st cluster
    {
        double u1 = -PI/2;
        double u2 = -PI/2;
        double s = 0.3;
        double period = 2*PI;
        for(int i=0; i < 4000; i++) {
            vector<double> point(2);
            point[0] = periodicGaussianSample(u1, s, period);
            point[1] = periodicGaussianSample(u2, s, period);
            dataset.push_back(point);
        }
    }

    // setup 2nd cluster
    {
        double u1 = PI/4;
        double u2 = 0;
        double s = 0.3;
        double period = 2*PI;
        for(int i=0; i < 4000; i++) {
            vector<double> point(2);
            point[0] = periodicGaussianSample(u1, s, period);
            point[1] = periodicGaussianSample(u2, s, period);
            dataset.push_back(point);
        }
    }

    // uncomment if visualization is desired
    // plot the output using xmgrace directly if you want to view it
    // ex. ./testAssignment > data.txt; xmgrace data.txt
    /*
    for(int i=0; i < dataset.size(); i++) {
        cout << dataset[i][0] << " " << dataset[i][1] << endl;
    }
    */

}

int main() {
    testEasyCase2D();

}
