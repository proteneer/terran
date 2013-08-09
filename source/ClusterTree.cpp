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
    period_(period),
    root_(NULL) {

    vector<int> points(dataset.size());
    for(int i=0; i<points.size(); i++) {
        points[i]=i;
    }
    root_ = new Node;
    root_->indices = points;
    queue_.push(root_);
}

ClusterTree::~ClusterTree() {
    delete root_;
}

int ClusterTree::getNumPoints() const {
    return dataset_.size();
}


/*
vector<int> ClusterTree::getPointsInCluster(int clusterIndex) const {
    return clusters_[clusterIndex].indices;
}

*/

/*
vector<int> ClusterTree::getAssignment() const {
    
    vector<int> allIndices(getNumPoints());
    int clusterIndex = 0;
    for(int i=0; i < clusters_.size(); i++) {
        vector<int> points = clusters_[i].indices;
        for(int j=0; j < points.size(); j++) {
            allIndices[points[i]] = clusterIndex;
        }
        clusterIndex++;
    }
    return allIndices;
    
}
*/

ofstream log2("log2.txt");

bool ClusterTree::step() {
    if(queue_.size() == 0) 
        return false;
    Node* node = queue_.front();
    queue_.pop();

    // nodes in the queue have the property:
    // 1. partitions are not set
    // 2. children do not exist
    if(node->indices.size() == 0) 
        throw(std::runtime_error("ClusterTree::step() - node has no points"));
    if(node->partitions.size() > 0)
        throw(std::runtime_error("ClusterTree::step() - node partition not empty"));
    if(node->children.size() > 0)
        throw(std::runtime_error("ClusterTree::step() - node children not empty"));
    
    // if this node doesn't have data points then break
    if(node->indices.size() < 500)
        return true;

    vector<vector<double> > subset;
    const vector<int> &indices = node->indices;
    for(int i=0; i<indices.size(); i++) {
        subset.push_back(dataset_[indices[i]]);
    }
    Cluster cluster(subset, period_);

    // setup partition tools

    vector<int> assignment = cluster.cluster();
    int maxAssignmentId = *(max_element(assignment.begin(), assignment.end()));

    if(maxAssignmentId > 0) {
        vector<vector<int> > subsetIndices(maxAssignmentId+1);
        for(int j=0; j < assignment.size(); j++) {
            subsetIndices[assignment[j]].push_back(indices[j]);
        }
        // for each new cluster
        for(int j=0; j < subsetIndices.size(); j++) {
            Node* newNode = new Node;
            newNode->indices = subsetIndices[j];
            // set the children for this cluster
            node->children.push_back(newNode);
        }
    }
    return true;
}

ClusterTree::Node* ClusterTree::getRoot() {
    return root_;
}