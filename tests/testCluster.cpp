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


void testClusterTree() {

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
    ClusterTree ct(dataset, periodset);

    while(!ct.finished()) {
        ct.setCurrentCluster();
        for(int i=0 ; i < ct.getCurrentCluster().getNumDimensions(); i++) {
            ct.partitionCurrentCluster(i);
            vector<double> partitions = ct.getCurrentCluster().getPartition(i);
        }
        ct.divideCurrentCluster();
    }

    vector<int> assignment = ct.getAssignment();
    ofstream l0("log0.txt");
    ofstream l1("log1.txt");
    ofstream l2("log2.txt");

    for(int i=0; i < assignment.size(); i++) {
        if(assignment[i] == 0) {
            l0 << ct.getPoint(i)[0] << " " << ct.getPoint(i)[1] << endl;
        }
        if(assignment[i] == 1) {
            l1 << ct.getPoint(i)[0] << " " << ct.getPoint(i)[1] << endl;
        }        
        if(assignment[i] == 2) {
            l2 << ct.getPoint(i)[0] << " " << ct.getPoint(i)[1] << endl;
        }

    }

    int numClusters = *(max_element(assignment.begin(), assignment.end()))+1;

    if(numClusters != 3) {
        throw(std::runtime_error("Wrong number of clusters!"));
    }

}

//  PI
//  |        **3** 
//  |*        ***       *        
//  |**                **
// d|0**              **0 
// i|**       ***      **
// m|*       **2**      * 
//  |        *****       
// 1|*        ***       *
//  |**                **
//  |1*       ***      *1
//  |*       **3**      *  
// -PI-------------------PI
//         dim 0

void testMultiCluster() {

    double period = 2*PI;
    vector<vector<double> > dataset;
    // setup 0th cluster;
    {
        double u1 = -4*PI/5;
        double u2 =  PI/2;
        double s1 = 0.4;
        double s2 = 0.7;
        
        for(int i=0; i < 2000; i++) {
            vector<double> point(2);
            point[0] = periodicGaussianSample(u1, s1, period);
            point[1] = periodicGaussianSample(u2, s2, period);
            dataset.push_back(point);
        }
    }

    // setup 1st cluster
    {
        double u1 = 7*PI/8;
        double u2 = -PI/2;
        double s1 = 0.4;
        double s2 = 0.6;

        for(int i=0; i < 2000; i++) {
            vector<double> point(2);
            point[0] = periodicGaussianSample(u1, s1, period);
            point[1] = periodicGaussianSample(u2, s2, period);
            dataset.push_back(point);
        }
    }

    // setup 2nd cluster
    {
        double u1 = PI/8;
        double u2 = -PI/9;
        double s = 0.5;

        for(int i=0; i < 2000; i++) {
            vector<double> point(2);
            point[0] = periodicGaussianSample(u1, s, period);
            point[1] = periodicGaussianSample(u2, s, period);
            dataset.push_back(point);
        }
    }

    // setup 3rd cluster
    {
        double u1 = 0;
        double u2 = 8*PI/9;
        double s1 = 0.4;
        double s2 = 0.7;

        for(int i=0; i < 2000; i++) {
            vector<double> point(2);
            point[0] = periodicGaussianSample(u1, s1, period);
            point[1] = periodicGaussianSample(u2, s2, period);
            dataset.push_back(point);
        }
    }

    vector<double> periodset;
    periodset.push_back(2*PI);
    periodset.push_back(2*PI);
    ClusterTree ct(dataset, periodset);



    while(!ct.finished()) {
        ct.setCurrentCluster();
        for(int i=0 ; i < ct.getCurrentCluster().getNumDimensions(); i++) {
            ct.partitionCurrentCluster(i);
            vector<double> partitions = ct.getCurrentCluster().getPartition(i);
            cout << "dim " << i << " ";
            for(int j=0; j < partitions.size(); j++) {
                cout << partitions[j] << " ";
            }
            cout << endl;
        }
        ct.divideCurrentCluster();
    }

/*
    ct.setCurrentCluster();
    for(int d=0 ; d < ct.getCurrentCluster().getNumDimensions(); d++) {
        ct.partitionCurrentCluster(d);
        vector<double> partitions = ct.getCurrentCluster().getPartition(d);
        cout << "dim " << d << " ";
        for(int i=0; i < partitions.size(); i++) {
            cout << partitions[i] << " ";
        }
        cout << endl;
    }
    ct.divideCurrentCluster();

    ct.setCurrentCluster();
    for(int d=0 ; d < ct.getCurrentCluster().getNumDimensions(); d++) {
        ct.partitionCurrentCluster(d);
        vector<double> partitions = ct.getCurrentCluster().getPartition(d);
        cout << "dim " << d << " ";
        for(int i=0; i < partitions.size(); i++) {
            cout << partitions[i] << " ";
        }
        cout << endl;
    }
    ct.divideCurrentCluster();
    */

    
    ofstream datalog("cdata0.txt");
    for(int i=0; i< dataset.size(); i++) {
        datalog << dataset[i][0] << " " << dataset[i][1] << endl;
    }
    

}



int main() {
    try{
        //testEasyCase2D();
        srand(1);
        //testClusterTree();
        srand(1);
        testMultiCluster();
        cout << "done" << endl;
    } catch(const exception &e) {
        cout << e.what();
    }
}
