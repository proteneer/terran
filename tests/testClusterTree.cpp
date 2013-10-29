#include <iostream>
#include <fstream>
#include <algorithm>

#include <Terran.h>
#include "util.h"

#include <sstream>

using namespace Terran;
using namespace std;


void testPeriodicSimpleCase() {
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
	double s2 = 0.4;
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
	double s2 = 0.4;
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
    ClusterTree ct(dataset, periodset);

	
    while(ct.queueSize() > 0) {
        ct.setCurrentCluster();
	ct.getCurrentCluster().partitionAll();
	ct.divideCurrentCluster(3000);
    }

    vector<int> assignment = ct.assign();

    int numClusters = *(max_element(assignment.begin(), assignment.end()))+1;

    if(numClusters != 3) {
        throw(std::runtime_error("testClusterTree() - Wrong number of clusters!"));
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
//
// Used in paper
void testPeriodicMultiCluster() {

    double period = 2*PI;
    vector<vector<double> > dataset;
    // setup 0th cluster;
    {
        double u1 = -4*PI/5;
        double u2 =  PI/2;
        double s1 = 0.3;
        double s2 = 0.45;
        
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
        double s1 = 0.3;
        double s2 = 0.45;

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
        double s = 0.4;

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
        double s1 = 0.3;
        double s2 = 0.5;

        for(int i=0; i < 2000; i++) {
            vector<double> point(2);
            point[0] = periodicGaussianSample(u1, s1, period);
            point[1] = periodicGaussianSample(u2, s2, period);
            dataset.push_back(point);
        }
    }

    vector<int> periodset(2,true);
    ClusterTree ct(dataset, periodset);

    while(ct.queueSize() > 0) {
        ct.setCurrentCluster();
	ct.getCurrentCluster().partitionAll();
	ct.divideCurrentCluster(3000);
    }

    vector<int> assignment = ct.assign();

    int numClusters = *(max_element(assignment.begin(), assignment.end()))+1;

    if(numClusters != 4) {
        throw(std::runtime_error("testPeriodicMultiCluster() - Wrong number of clusters!"));
    }

}

int main() {
    try{
        cout << "testPeriodicSimpleCase()" << endl;
	srand(1);
        testPeriodicSimpleCase();
        srand(1);
        cout << "testPeriodicMultiCluster()" << endl;
        testPeriodicMultiCluster();
        //cout << "testNonPeriodicMultiCluster()" << endl;
        //srand(1);
        //testNonPeriodicMultiCluster();
        cout << "done" << endl;
    } catch(const exception &e) {
        cout << e.what();
    }
}
