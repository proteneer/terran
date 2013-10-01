#include <iostream>
#include <fstream>
#include <algorithm>

#include <Terran.h>
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
void testPeriodicSimpleCase() {

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

    while(!ct.finished()) {
        ct.setCurrentCluster();
        for(int i=0 ; i < ct.getCurrentCluster().getNumDimensions(); i++) {
            ct.partitionCurrentCluster(i);
            vector<double> partitions = ct.getCurrentCluster().getPartition(i);
        }

        // if we use a non reference, and we do something like
        // Cluster cc = ct.getCurrentCluster();
        // we will end up with a segfault as the references members get destroyed
        Cluster& cc = ct.getCurrentCluster();
        vector<vector<double> > clusterPoints;
        for(int i=0 ; i < cc.getNumPoints(); i++) {
            clusterPoints.push_back(cc.getPoint(i));
        }
        vector<int> assign = cc.cluster();
        ct.divideCurrentCluster();
        count++;
    }

    vector<int> assignment = ct.getAssignment();

    int numClusters = *(max_element(assignment.begin(), assignment.end()))+1;

    if(numClusters != 4) {
        throw(std::runtime_error("testPeriodicMultiCluster() - Wrong number of clusters!"));
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
// used for the paper
/*
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

        Cluster& cc = ct.getCurrentCluster();

        for(int i=0 ; i < ct.getCurrentCluster().getNumDimensions(); i++) {
            ct.partitionCurrentCluster(i);

            cout << "dim " << i << " ";

            // PLOTTING
            cout << "gaussian components" << endl;
            try {

                PartitionerEM &pem = dynamic_cast<PartitionerEM&>(cc.getPartitioner(i));
                vector<Param> params = pem.getEM().getParams();

                stringstream filename;
                filename << "curves" << count << i;

                if(i == 0) {
                    Util::plotGaussian(params, -2, 4.5, filename.str());
                } else if (i == 1) {
                    Util::plotGaussian(params, -3.9, 5, filename.str());
                }

            } catch(...) {} 
            vector<double> partitions = ct.getCurrentCluster().getPartition(i);


            for(int j=0; j < partitions.size(); j++) {
                cout << partitions[j] << " ";
            }
            cout << endl;
        }

        // if we use a non reference, and we do something like
        // Cluster cc = ct.getCurrentCluster();
        // we will end up with a segfault as the references get destroyed


        vector<vector<double> > clusterPoints;
        for(int i=0 ; i < cc.getNumPoints(); i++) {
            clusterPoints.push_back(cc.getPoint(i));
        }
        vector<int> assign = cc.cluster();

        stringstream filename;
        filename << "X" << count;
        ofstream X(filename.str().c_str());

        for(int i=0; i < cc.getNumPoints(); i++) {
 
            X << clusterPoints[i][0] << " " << clusterPoints[i][1] << endl;


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

    ofstream assignmentFile("assignment.log");
    for(int i=0; i< assignment.size() ;i++) {
        assignmentFile << assignment[i] << endl;
    }
}
*/

int main() {
    try{
        cout << "testPeriodicSimpleCase()" << endl;
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
