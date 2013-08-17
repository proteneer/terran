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

using namespace Terran;
using namespace std;

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
void testPeriodicMultiCluster() {

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

    int count = 0;

    //while(count < 2) {
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

        // if we use a non reference, and we do something like
        // Cluster cc = ct.getCurrentCluster();
        // we will end up with a segfault as the references get destroyed

        Cluster& cc = ct.getCurrentCluster();
        vector<vector<double> > clusterPoints;
        for(int i=0 ; i < cc.getNumPoints(); i++) {
            clusterPoints.push_back(cc.getPoint(i));
        }
        vector<int> assign = cc.cluster();

        for(int i=0; i < cc.getNumPoints(); i++) {
            stringstream filename;
            int j = assign[i];
            filename << "paper" << count << j;
            ofstream fname(filename.str().c_str(), ios::app);
            fname << clusterPoints[i][0] << " " << clusterPoints[i][1] << endl;
            fname.close();
        }
        ct.divideCurrentCluster();
        count++;
    }

    vector<int> assignment = ct.getAssignment();

    for(int i=0; i < ct.getNumPoints(); i++) {
        stringstream filename;
        int j = assignment[i];
        filename << "final" << count << j;
        ofstream fname(filename.str().c_str(), ios::app);
        fname << dataset[i][0] << " " << dataset[i][1] << endl;
        fname.close();
    }

    ofstream datalog("cdata0.txt");
    for(int i=0; i< dataset.size(); i++) {
        datalog << dataset[i][0] << " " << dataset[i][1] << endl;
    }

}
//+inf
//  |         ***
//  |        **3** 
//  |         ***       *        
//  |                  ***
// d|                 **0** 
// i|         ***      ***
// m|        **2**      * 
//  |        *****       
// 1|         ***       *
//  |                  ***
//  |                  *1*
//  |                   *  
// -PI-------------------PI--+inf
//  |
//-inf       dim 0
void testNonPeriodicMultiCluster() {

    vector<vector<double> > dataset;
    // setup 0th cluster;
    {
        double u1 =  4*PI/5;
        double u2 =  PI/2;
        double s1 = 0.4;
        double s2 = 0.7;
        
        for(int i=0; i < 2000; i++) {
            vector<double> point(2);
            point[0] = gaussianSample(u1, s1);
            point[1] = gaussianSample(u2, s2);
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
            point[0] = gaussianSample(u1, s1);
            point[1] = gaussianSample(u2, s2);
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
            point[0] = gaussianSample(u1, s);
            point[1] = gaussianSample(u2, s);
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
            point[0] = gaussianSample(u1, s1);
            point[1] = gaussianSample(u2, s2);
            dataset.push_back(point);
        }
    }


    ofstream datalog("cdata0.txt");
    for(int i=0; i< dataset.size(); i++) {
        datalog << dataset[i][0] << " " << dataset[i][1] << endl;
    }

    vector<double> periodset;
    periodset.push_back(0);
    periodset.push_back(0);
    ClusterTree ct(dataset, periodset);

    int count = 0;

    //while(count < 2) {
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

        // if we use a non reference, and we do something like
        // Cluster cc = ct.getCurrentCluster();
        // we will end up with a segfault as the references get destroyed

        Cluster& cc = ct.getCurrentCluster();
        vector<vector<double> > clusterPoints;
        for(int i=0 ; i < cc.getNumPoints(); i++) {
            clusterPoints.push_back(cc.getPoint(i));
        }
        vector<int> assign = cc.cluster();

        for(int i=0; i < cc.getNumPoints(); i++) {
            stringstream filename;
            int j = assign[i];
            filename << "paper" << count << j;
            ofstream fname(filename.str().c_str(), ios::app);
            fname << clusterPoints[i][0] << " " << clusterPoints[i][1] << endl;
            fname.close();
        }
        ct.divideCurrentCluster();
        count++;
    }

    vector<int> assignment = ct.getAssignment();

    for(int i=0; i < ct.getNumPoints(); i++) {
        stringstream filename;
        int j = assignment[i];
        filename << "final" << count << j;
        ofstream fname(filename.str().c_str(), ios::app);
        fname << dataset[i][0] << " " << dataset[i][1] << endl;
        fname.close();
    }


}

int main() {
    try{
        //testClusterTree();
        //srand(1);
        //testMultiCluster();
        srand(1);
        testNonPeriodicMultiCluster();
        cout << "done" << endl;
    } catch(const exception &e) {
        cout << e.what();
    }
}
