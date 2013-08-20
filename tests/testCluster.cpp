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


void test100DGaussians() {

    const int numDimensions = 100;
    const int numModes = 3;

    // generate 100 clusters with different means
    vector<vector<double> > dataset;
    vector<double> us(numModes);
    vector<double> ss(numModes);

    us[0] = 0;
    ss[0] = 0.35;

    us[1] = 4;
    ss[1] = 0.7;

    us[2] = 8;
    ss[2] = 0.5;

/*
    for(int d=0; d < numModes; d++) {
        us[d] = d*2;
        ss[d] = 0.35;
    }
*/
    // initialize a set of clusters
    /*
    for(int d=0; d<numModes; d++) {
        for(int i=0; i < 1000; i++) {
            vector<double> point(numDimensions);
            for(int j = 0; j < numDimensions; j++) {
                point[j] = gaussianSample(us[numModes], ss[j%numModes]);
            }
            dataset.push_back(point);
        }
    }
    */

    for(int k=0; k < 10; k++) {
        for(int i=0; i < 1000; i++) {
            vector<double> point(numDimensions);
            for(int j = 0; j < numDimensions; j++) {
                if( j == k )
                    point[j] = gaussianSample(us[0], ss[0]);
                else
                    point[j] = gaussianSample(us[j%numModes], ss[j%numModes]);
            }
            dataset.push_back(point);
        }   
    }

    for(int i=0; i < 1000; i++) {
        vector<double> point(numDimensions);
        for(int j = 0; j < numDimensions; j++) {
            point[j] = gaussianSample(us[j%numModes], ss[j%numModes]);
        }
        dataset.push_back(point);
    }

/*
    for(int i=0; i < 1000; i++) {
        vector<double> point(numDimensions);
        for(int j = 0; j < numDimensions; j++) {
            if( j == 0 )
                point[j] = gaussianSample(us[2], ss[2]);
            else
                point[j] = gaussianSample(us[j], ss[j]);
        }
        dataset.push_back(point);
    }
    */

    /*
    for(int k=0; k < 5; k++) {
        vector<double> udd(numDimensions);
        vector<double> sdd(numDimensions);

        for(int d = 0; d < numDimensions; d++) {
            int sample = rand()%numDimensions;
            udd[d] = us[sample];
            sdd[d] = ss[sample];
        }
        for(int i=0; i < 1000; i++) {
            // decide where to draw each dimension from
            
            vector<double> point(numDimensions);

            for(int d=0; d < numDimensions; d++) {
                point[d] = gaussianSample(udd[d], sdd[d]);
            }
            dataset.push_back(point);
        }
    }
    */

    /*
    for(int j=0; j < 3; j++) {
        vector<double> point(numDimensions);
        for(int j = 0; j < numDimensions; j++) {
            double u = us[rand()%numDimensions];
            double s = us[rand()%numDimensions];
            for(int i=0; i < 1000; i++) {
                point[j] = gaussianSample(u, s);
            }
        }
        dataset.push_back(point);
    }*/


    // randomize the points
    random_shuffle(dataset.begin(), dataset.end());

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


    vector<double> d0 = cc.getDimension(0);

    ofstream ld0("d0.txt");

    for(int i=0; i <d0.size(); i++) {
        ld0 << d0[i] << endl;
    }




    vector<vector<double> > testPartitions;
    for(int d = 0; d < cc.getNumDimensions(); d++) {
        cout << "partitioning " << d << endl;
        cc.partition(d);
        vector<double> partitions = cc.getPartition(d);
        testPartitions.push_back(partitions);
    }

    vector<int> assignment = cc.cluster();

    ofstream assignmentd("100dassign.txt");
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



    if(numClusters != 3) {
        throw(std::runtime_error("test50DGaussians(): Number of clusters != 3"));
    }

}

/*
// generate three 50 dimensional gaussians
void test100DGaussians() {
    const int numDimensions = 100;
    vector<vector<double> > dataset;
    {
        double u = 3;
        double s = 0.3;
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

    {
        double u1 = 3;
        double s1 = 0.3;
        double u2 = -6;
        double s2 = 0.7;
        double u3 = 12;
        double s3 = 1.2;
        for(int i=0; i < 2000; i++) {
            vector<double> point(numDimensions);
            for(int d=0; d < numDimensions; d++) {
                if(d % 3 == 0)
                    point[d] = gaussianSample(u1,s1);
                else if(d % 3 == 1)
                    point[d] = gaussianSample(u2,s2);
                else
                    point[d] = gaussianSample(u3,s3);
            }
            dataset.push_back(point);
        }
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
        cc.partition(d);
        vector<double> partitions = cc.getPartition(d);
        testPartitions.push_back(partitions);
    }

    vector<int> assignment = cc.cluster();

    ofstream assignmentd("100dassign.txt");
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
*/

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
        test100DGaussians();
        srand(1);
        testEasyCase2D();
        cout << "done" << endl;
    } catch(const exception &e) {
        cout << e.what();
    }
}
