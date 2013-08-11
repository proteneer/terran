#include <math.h>
#include <vector>
#include <iostream>
#include <fstream>

#include <EMGaussian.h>
#include <EMPeriodicGaussian.h>
#include <MethodsPeriodicGaussian.h>
#include <MathFunctions.h>
#include <Cluster.h>

#include <ClusterTree.h>

#include "util.h"


#include <sstream>

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

using namespace Terran;


void testEasyCase2D() {

    vector<vector<double> > dataset;
    // setup 0th cluster;
    {
        double u1 = -PI/2;
        double u2 =  PI/2;
        double s = 0.5;
        double period = 2*PI;
        for(int i=0; i < 2000; i++) {
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
        double s = 0.5;
        double period = 2*PI;
        for(int i=0; i < 2000; i++) {
            vector<double> point(2);
            point[0] = periodicGaussianSample(u1, s, period);
            point[1] = periodicGaussianSample(u2, s, period);
            dataset.push_back(point);
        }
    }

    // setup 2nd cluster
    {
        double u1 = PI/3;
        double u2 = 0;
        double s = 0.5;
        double period = 2*PI;
        for(int i=0; i < 2000; i++) {
            vector<double> point(2);
            point[0] = periodicGaussianSample(u1, s, period);
            point[1] = periodicGaussianSample(u2, s, period);
            dataset.push_back(point);
        }
    }


    vector<double> periodset;
    periodset.push_back(2*PI);
    periodset.push_back(2*PI);

    Cluster cc(dataset, periodset);

    ofstream fdata("fdata.txt");
    for(int i=0; i < dataset.size(); i++) {
        fdata << dataset[i][0] << " " << dataset[i][1] << endl;
    }

    vector<vector<double> > testPartitions;
    
    for(int d = 0; d < cc.getNumDimensions(); d++) {
        cc.partition(d);
        vector<double> partitions = cc.getPartition(d);
        testPartitions.push_back(partitions);
    }

    vector<double> truthPartition0;
    vector<double> truthPartition0Errors;
    truthPartition0.push_back(-0.10);
    truthPartition0.push_back( 2.82);
    truthPartition0Errors.push_back(0.1);
    truthPartition0Errors.push_back(0.2);

    vector<double> truthPartition1;
    vector<double> truthPartition1Errors;
    // +- 0.15 is near the periodic boundary, so we need to test with matchPeriodicPoints
    truthPartition1.push_back(-3.1);
    truthPartition1Errors.push_back(0.15);

    for(int i=0; i < testPartitions.size(); i++) {
        for(int j=0; j < testPartitions[i].size(); j++) {
            cout << testPartitions[i][j] << " ";
        }
        cout << endl;
    }

    /*
    ofstream cluster0("cdata0.txt");
    ofstream cluster1("cdata1.txt");

    vector<int> assignment = cc.cluster();
    for(int i=0; i<assignment.size(); i++) {
        if(assignment[i] == 0) {
            cluster0 << dataset[i][0] << " " << dataset[i][1] << endl;
        } else if(assignment[i] == 1) {
            cluster1 << dataset[i][0] << " " << dataset[i][1] << endl;
        }
    }
    */

    Util::matchPoints(truthPartition0, testPartitions[0], truthPartition0Errors);
    Util::matchPeriodicPoints(truthPartition1, testPartitions[1], 2*PI, truthPartition0Errors);

}

/*
int testClusterTree() {

    vector<vector<double> > dataset;
    vector<double> periodset;
    // setup 0th cluster;
    {
        double u1 = -PI/2;
        double u2 =  PI/2;
        double s = 0.5;
        double period = 2*PI;
        for(int i=0; i < 2000; i++) {
            vector<double> point(2);
            point[0] = periodicGaussianSample(u1, s, period);
            point[1] = periodicGaussianSample(u2, s, period);
            dataset.push_back(point);
            periodset.push_back(period);
        }
    }

    // setup 1st cluster
    {
        double u1 = -PI/2;
        double u2 = -PI/2;
        double s = 0.5;
        double period = 2*PI;
        for(int i=0; i < 2000; i++) {
            vector<double> point(2);
            point[0] = periodicGaussianSample(u1, s, period);
            point[1] = periodicGaussianSample(u2, s, period);
            dataset.push_back(point);
            periodset.push_back(period);
        }
    }

    // setup 2nd cluster
    {
        double u1 = PI/3;
        double u2 = 0;
        double s = 0.5;
        double period = 2*PI;
        for(int i=0; i < 2000; i++) {
            vector<double> point(2);
            point[0] = periodicGaussianSample(u1, s, period);
            point[1] = periodicGaussianSample(u2, s, period);
            dataset.push_back(point);
            periodset.push_back(period);
        }
    }

    ClusterTree ct(dataset, periodset);
    ct.stepBFS();
    ct.stepBFS();
    ct.stepBFS();
}
*/

int main() {
    try{
        testEasyCase2D();
        //testClusterTree();
        cout << "done" << endl;
    } catch(const exception &e) {
        cout << e.what();
    }
}
