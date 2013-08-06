#include "ClusterTree.h"
#include "Cluster.h"
#include <map>
#include <iostream>
#include <utility>

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

}

ClusterTree::~ClusterTree() {

}

// go down one level
bool ClusterTree::stepBFS() {

    // TODO: make sure clusters_ on each loop have every single point

    // the size of clusters_ will change during the loop, so we keep track of how many trips to make. 
    // new clusters are appended to the end of clusters_, so it is OK to overwrite them.
    int initialSize = clusters_.size();

    for(int i=0; i < initialSize; i++) {
        if(clusters_[i].second == false) {
            vector<int> clusterPointIndices = clusters_[i].first;
            // generate the point data for this cluster
            vector<vector<double> > clusterPointData;
            for(int j=0; j < clusterPointIndices.size(); j++) {
                clusterPointData.push_back(dataset_[clusterPointIndices[j]]);
            }

            // optimize the parameters for each dimension
            Cluster cc(clusterPointData, period_);
            for(int d=0; d < cc.getNumDimensions(); d++) {
                cc.optimizeParameters(d);
                cc.partition(d,0.05);
            }
            // ranged from 0-(C-1), where C is the number of clusters
            vector<int> assignment = cc.cluster();

            // determine the max number of clusters
            int maxAssignmentId = 0;
            for(int j=0; j < assignment.size(); j++) {
                maxAssignmentId = assignment[j] > maxAssignmentId ? assignment[j] : maxAssignmentId;
            }

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
}


// go down one level
/*
void ClusterTree::BFS() {
    
    for(int i=0; i < clusters_.size(); i++) {
            vector<int> clusterPointIndices = clusters_[i].first;
            // generate the point data for this cluster
            vector<vector<double> > clusterPointData;
            for(int j=0; j < clusterPointIndices.size(); j++) {
                clusterPointData.push_back(dataset_[clusterPointIndices[j]]);
            }

            // optimize the parameters for each dimension
            Cluster cc(clusterPointData, period_);
            for(int d=0; d < cc.getNumDimensions(); d++) {
                cc.optimizeParameters(d);
                cc.partition(d,0.05);
            }
            // ranged from 0-(C-1), where C is the number of clusters
            vector<int> assignment = cc.cluster();

            // determine the max number of clusters
            int maxAssignmentId = 0;
            for(int j=0; j < assignment.size(); j++) {
                maxAssignmentId = assignment[j] > maxAssignmentId ? assignment[j] : maxAssignmentId;
            }
            clusters_[i].first = newClusterPointIndices[0];
            // if there exists more than one 

            if(maxAssignmentId > 0) {

                vector<vector<int> > newClusterPointIndices(maxAssignmentId+1);
                for(int j=0; j < assignment.size(); j++) 
                    newClusterPointIndices[assignment[j]].push_back(clusterPointIndices[j]);
            
                // partition up the data
                
                for(int j=1; j < newClusterPointIndices.size(); j++)
                    clusters_.push_back(make_pair(newClusterPointIndices[j],false));

            // otherwise this cluster cannot be partitioned any further
            } else {
                clusters_[i].second = true;
            }
        }
}
*/