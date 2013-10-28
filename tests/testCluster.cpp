#include <math.h>
#include <vector>
#include <iostream>
#include <fstream>
#include <algorithm>

#include <Terran.h>

#include "util.h"

#include <sstream>

using namespace std;
using namespace Terran;

void testEasyCase2D() {

	//  three-cluster periodic dataset:

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

    vector<vector<double> > dataset;
    // setup 0th cluster;
    {
        double u1 = -PI/2;
        double u2 =  PI/2;
        double s1 = 0.3;
		double s2 = 0.5;
        double period = 2*PI;
        for(int i=0; i < 2000; i++) {
            vector<double> point(2);
            point[0] = periodicGaussianSample(u1, s1, period);
            point[1] = periodicGaussianSample(u2, s2, period);
            dataset.push_back(point);
        }
    }

    // setup 1st cluster
    {
        double u1 = -PI/2;
        double u2 = -PI/2;
        double s1 = 0.3;
		double s2 = 0.5;
        double period = 2*PI;
        for(int i=0; i < 2000; i++) {
            vector<double> point(2);
            point[0] = periodicGaussianSample(u1, s1, period);
            point[1] = periodicGaussianSample(u2, s2, period);
            dataset.push_back(point);
        }
    }

    // setup 2nd cluster
    {
        double u1 = PI/2;
        double u2 = 0;
        double s1 = 0.3;
		double s2 = 0.5;
        double period = 2*PI;
        for(int i=0; i < 2000; i++) {
            vector<double> point(2);
            point[0] = periodicGaussianSample(u1, s1, period);
            point[1] = periodicGaussianSample(u2, s2, period);
            dataset.push_back(point);
        }
    }

    vector<int> periodset(2, true);

    Cluster cc(dataset, periodset);

    vector<vector<double> > testPartitions;
    
    for(int d = 0; d < cc.getNumDimensions(); d++) {
        cc.partition(d);
        vector<double> partitions = cc.getPartition(d);
        testPartitions.push_back(partitions);
    }

    vector<double> truthPartition0;
    vector<double> truthPartition0Errors;
    truthPartition0.push_back(0.0);
    truthPartition0.push_back(3.0);
    truthPartition0Errors.push_back(0.1);
    truthPartition0Errors.push_back(0.2);

    vector<double> truthPartition1;
    vector<double> truthPartition1Errors;

    truthPartition1.push_back(-3.14);
    truthPartition1Errors.push_back(0.1);

    Util::matchPoints(truthPartition0, testPartitions[0], truthPartition0Errors);
    Util::matchPeriodicPoints(truthPartition1, testPartitions[1], 2*PI, truthPartition1Errors);
}

int main() {
    try{
        testEasyCase2D();
        cout << "done" << endl;
    } catch(const exception &e) {
        cout << e.what();
    }
}
