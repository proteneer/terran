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

/*
void test100DGaussiansMulti() {

    const int numDimensions = 100;
    const int numModes = 6;

    // generate 100 clusters with different means
    vector<vector<double> > dataset;
    vector<double> us(numModes);
    vector<double> ss(numModes);

    // generate randomly
    for(int i=0; i < numModes; i++) {
        us[i] = 2*i;
        if(i%2 == 0)
            ss[i] = 0.25;
        else
            ss[i] = 0.15;
    }

    ofstream flabels_("100dMultiLabels.txt");

    for(int k=0; k < 20; k++) {

        vector<double> means(numDimensions);
        vector<double> stddev(numDimensions);

        for(int d=0; d < numDimensions; d++) {
            int sample = rand() % numModes;
            means[d] = us[sample];
            stddev[d] = ss[sample];
        }

        for(int i=0; i < 1000; i++) {
            vector<double> point(numDimensions);
            for(int d = 0; d < numDimensions; d++) {
                point[d] = gaussianSample(means[d], stddev[d]);
            }
            dataset.push_back(point);
            flabels_ << k << endl;
        }
    }

    // k means
    ofstream highd("100dMulti.txt");
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
        cout << "partitioning " << d << "...";
        cc.partition(d);
        vector<double> partitions = cc.getPartition(d);
        cout << " found: " << partitions.size() << endl;
        testPartitions.push_back(partitions);
    }

    vector<int> assignment = cc.cluster();

    ofstream assignmentd("100dMultiTerranAssign.txt");
    for(int j=0; j < assignment.size(); j++) {
        assignmentd << assignment[j] << endl;
    }

    int numClusters = (*max_element(assignment.begin(), assignment.end()))+1;

    cout << "number of clusters: " << numClusters << endl;

    vector<vector<int> > output(numClusters);
    for(int i=0; i < assignment.size(); i++) {
        output[assignment[i]].push_back(i);
    }

    int realc = 0;
    for(int i=0; i < output.size(); i++) {
        if(output[i].size() > 100) {
            cout << output[i].size() << endl;
            realc++;
        }
    }

    cout << "number of real clusters: " << realc << endl;


}

void test100DGaussians() {

    const int numDimensions = 100;
    const int numModes = 3;

    // generate 100 clusters with different means
    vector<vector<double> > dataset;
    vector<double> us(numModes);
    vector<double> ss(numModes);
    
    ofstream flabels_("100dLabels.txt");

    // generate first set of points

    for(int i=0; i < 1000; i++) {
        // u= [0,0,...,0]
        vector<double> point(numDimensions);
        for(int d=0; d < numDimensions; d++) {
            if(d == 0) {
                point[d] = gaussianSample(0, 0.2);
            } else {
                point[d] = gaussianSample(0, 1.0);

            }
        }
        dataset.push_back(point);
        flabels_ << 0 << endl;
    }

    for(int i=0; i < 1000; i++) {
        // su = [-1,-1,...,-1]
        vector<double> point(numDimensions);
        for(int d=0; d < numDimensions; d++) {
            point[d] = gaussianSample(-1, 0.15);
        }
        dataset.push_back(point);
        flabels_ << 1 << endl;
    }

    ofstream highd("100d.txt");
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
        cout << "partitioning " << d << "...";
        cc.partition(d);
        vector<double> partitions = cc.getPartition(d);
        cout << " found: " << partitions.size() << endl;
        testPartitions.push_back(partitions);
    }

    vector<int> assignment = cc.cluster();

    ofstream assignmentd("100dTerranAssign.txt");
    for(int j=0; j < assignment.size(); j++) {
        assignmentd << assignment[j] << endl;
    }

    int numClusters = (*max_element(assignment.begin(), assignment.end()))+1;

    cout << "number of clusters: " << numClusters << endl;

    vector<vector<int> > output(numClusters);
    for(int i=0; i < assignment.size(); i++) {
        output[assignment[i]].push_back(i);
    }


    int realc = 0;
    for(int i=0; i < output.size(); i++) {
        if(output[i].size() > 100) {
            cout << output[i].size() << endl;
            realc++;
        }
    }

    cout << "number of real clusters: " << realc << endl;


}
*/

// A three-cluster periodic dataset looking something like this:

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

    vector<bool> periodset(2, true);

    Cluster cc(dataset, periodset);

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
    truthPartition0Errors.push_back(0.15);
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
        //test100DGaussiansMulti();
        //srand(1);
        //test100DGaussians();
        srand(1);
        testEasyCase2D();
        cout << "done" << endl;
    } catch(const exception &e) {
        cout << e.what();
    }
}
