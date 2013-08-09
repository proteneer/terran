#include "ClusterTree.h"
#include "Cluster.h"
#include <iostream>
#include <utility>
#include <algorithm>
#include <fstream>
#include <stdexcept>

using namespace std;
using namespace Terran;

ClusterTree::ClusterTree(const vector<vector<double> > &dataset, const vector<double> &period) : 
    dataset_(dataset),
    period_(period) {
    
    vector<int> points(dataset.size());
    for(int i=0; i<points.size(); i++) {
        points[i]=i;
    }
    vector<pair<vector<int>,bool> > initialClusters;
    initialClusters.push_back(make_pair(points,false));
    clusters_ = initialClusters;
}

ClusterTree::~ClusterTree() {

}

int ClusterTree::getNumPoints() const {
    return dataset_.size();
}

int ClusterTree::getNumClusters() const {
    return clusters_.size();
}

vector<int> ClusterTree::getPointsInCluster(int clusterIndex) const {
    return clusters_[clusterIndex].first;
}

vector<int> ClusterTree::getAssignment() const {
    vector<int> allIndices(getNumPoints());
    int clusterIndex = 0;
    for(int i=0; i < clusters_.size(); i++) {
        vector<int> points = clusters_[i].first;
        for(int j=0; j < points.size(); j++) {
            allIndices[points[i]] = clusterIndex;
        }
        clusterIndex++;
    }
    return allIndices;
}

ofstream log2("log2.txt");

// go down one level
bool ClusterTree::stepBFS() {

    // make sure clusters_ on each loop have every single point
    vector<int> allIndices;
    for(int i=0; i < clusters_.size(); i++) {
        allIndices.insert(allIndices.end(),(clusters_[i].first).begin(), (clusters_[i].first).end());
    }
    sort(allIndices.begin(), allIndices.end());
    if(allIndices.size() != getNumPoints()) {
        for(int i=0; i < allIndices.size(); i++) {
            cout << allIndices[i] << " ";
        }
        cout << allIndices.size() << " " << getNumPoints() << endl;
        throw(std::runtime_error("ClusterTree corruption: wrong number of points in clusters_"));
    }
    for(int i=0; i < allIndices.size()-1; i++) {
        if(allIndices[i] != allIndices[i+1] - 1) {
            throw(std::runtime_error("ClusterTree corruption: missing points in clusters_"));   
        } 
    }

    // the size of clusters_ will change during the loop, so we keep track of how many trips to make. 
    // new clusters are appended to the end of clusters_, so it is OK to overwrite them.
    int initialSize = clusters_.size();

    bool complete = true;



    log2 << "Clusters and sizes" << endl;
    for(int i=0; i < initialSize; i++) {
        log2 << clusters_[i].first.size() << " " << clusters_[i].second << endl;
    }

    for(int i=0; i < initialSize; i++) {
        log2 << "Working on cluster " << i << endl;
        if(clusters_[i].second == false) {
            complete = false;
            vector<int> clusterPointIndices = clusters_[i].first;
            // generate the point data for this cluster
            vector<vector<double> > clusterPointData;
            for(int j=0; j < clusterPointIndices.size(); j++) {
                clusterPointData.push_back(dataset_[clusterPointIndices[j]]);
            }

            // optimize the parameters for each dimension

            // alternative API:
            // for each dimension
            // 1. data = Cluster::getDimension(d); // get marginalized data for the d'th dimension
            // 2. vec2 = Partition(data); (either EM Partition or GKDE Partition)
            // 3. Cluster::setPartition(d);
            // 4. Cluster::cluster(); 

            // this way the Cluster class doesn't have to know about parameters, only the partitions;
            Cluster cc(clusterPointData, period_);

            /*
            // TODO FIGURE OUT HOW TO SET PARAMETERS HERE!
            for(int d=0; d < cc.getNumDimensions(); d++) {
                cc.optimizeParameters(d);
                cc.partition(d,0.05);
            }
            */

            log2 << "partitions: " << endl;
            vector<double> p0 = cc.getPartitions(0);
            for(int j=0; j<p0.size(); j++) {
                log2 << p0[j] << " ";
            }
            log2 << endl;

            vector<double> p1 = cc.getPartitions(1);
            for(int j=0; j<p1.size(); j++) {
                log2 << p1[j] << " ";
            }
            log2 << endl;

            // ranged from 0-(C-1), where C is the number of clusters
            vector<int> assignment = cc.cluster();

            // determine the max number of clusters
            int maxAssignmentId = 0;

            cout << assignment.size() << endl;

            for(int j=0; j < assignment.size(); j++) {
                log2 << assignment[j] << "";
                maxAssignmentId = assignment[j] > maxAssignmentId ? assignment[j] : maxAssignmentId;
            }
            log2 << endl;

            // if there exists more than one cluster
            if(maxAssignmentId > 0) {
                vector<vector<int> > newClusterPointIndices(maxAssignmentId+1);
                for(int j=0; j < assignment.size(); j++) {
                    newClusterPointIndices[assignment[j]].push_back(clusterPointIndices[j]);
                }
                for(int j=0; j < newClusterPointIndices.size(); j++) {
                    if(j == 0) {
                        // overwrite this cluster
                        clusters_[i] = make_pair(newClusterPointIndices[j],false);
                    } else {
                        // append to the end
                        clusters_.push_back(make_pair(newClusterPointIndices[j],false));
                    }
                }
            // otherwise mark as true
            } else {
                clusters_[i].second = true;
            }
        }
    }

    log2 << "Clusters and sizes" << endl;
    for(int i=0; i < clusters_.size(); i++) {
        log2 << clusters_[i].first.size() << " " << clusters_[i].second << endl;
    }

    return complete;
}
