#include <math.h>
#include <vector>
#include <iostream>
#include <fstream>
#include <algorithm>

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

// generate three 50 dimensional gaussians
void test50DGaussians() {
    const int numDimensions = 100;
    vector<vector<double> > dataset;
    {
        double u = 3;
        double s = 0.2;
        for(int i=0; i < 2000; i++) {
            vector<double> point(numDimensions);
            for(int d=0; d < numDimensions; d++) {
                point[d] = gaussianSample(u,s);
            }
            dataset.push_back(point);
        }
    }

    {
        double u = -6;
        double s = 0.7;
        for(int i=0; i < 2000; i++) {
            vector<double> point(numDimensions);
            for(int d=0; d < numDimensions; d++) {
                point[d] = gaussianSample(u,s);
            }
            dataset.push_back(point);
        }
    }

    // flip flop
    {
        double u = 12;
        double s = 1.2;
        for(int i=0; i < 2000; i++) {
            vector<double> point(numDimensions);
            for(int d=0; d < numDimensions; d++) {
                point[d] = gaussianSample(u,s);
            }
            dataset.push_back(point);
        }
    }

    ofstream highd("50d.txt");
    for(int i=0; i < dataset.size(); i++) {
        for(int j=0; j < dataset[i].size(); j++) {
            highd << dataset[i][j];
            if( j != dataset[i].size()-1) {
                highd << " ";
            }
        }
        highd << endl;
    }

    vector<double> periodset(numDimensions, 0);
    Cluster cc(dataset, periodset);
    vector<vector<double> > testPartitions;
    for(int d = 0; d < cc.getNumDimensions(); d++) {
        cc.partition(d);
        vector<double> partitions = cc.getPartition(d);
        testPartitions.push_back(partitions);
    }

    vector<int> assignment = cc.cluster();

    ofstream assignmentd("50dassign.txt");
    for(int j=0; j < assignment.size(); j++) {
        assignmentd << assignment[j] << endl;
    }

    int numClusters = (*max_element(assignment.begin(), assignment.end()))+1;

    cout << numClusters << endl;

    vector<vector<int> > output(numClusters);
    for(int i=0; i < assignment.size(); i++) {
        output[assignment[i]].push_back(i);
    }


    for(int i=0; i < output.size(); i++) {
        cout << output[i].size() << endl;
    }



    if(numClusters != 3) {
        throw(std::runtime_error("test50DGaussians(): Number of clusters != 3"));
    }

}

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

    vector<vector<double> > testPartitions;
    
    for(int d = 0; d < cc.getNumDimensions(); d++) {
        PartitionerEM& emp = dynamic_cast<PartitionerEM &>(cc.getPartitioner(d));
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

int main() {
    try{
        test50DGaussians();
        srand(1);
        testEasyCase2D();
        cout << "done" << endl;
    } catch(const exception &e) {
        cout << e.what();
    }
}
