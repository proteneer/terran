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
    root_(NULL),
    currentCluster_(NULL) {

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
    delete currentCluster_;
}

int ClusterTree::getNumPoints() const {
    return dataset_.size();
}

bool ClusterTree::finished() const {
    if(queue_.size() == 0) {
        return true;
    } else {
        return false;
    }
}


void ClusterTree::setCurrentCluster() {

    if(queue_.size() == 0) 
        throw(std::runtime_error("ClusterTree::step - invoked setCluster() on an empty queue"));
    if(currentCluster_ != NULL)
        throw(std::runtime_error("ClusterTree::currentCluster_ is not set to NULL, has divideCluster() been called?"));
    currentNode_ = queue_.front();
    queue_.pop();

    if(currentNode_->indices.size() == 0) 
        throw(std::runtime_error("ClusterTree::step() - currentNode_ has no points"));
    if(currentNode_->partitions.size() > 0)
        throw(std::runtime_error("ClusterTree::step() - currentNode_ partition not empty"));
    if(currentNode_->children.size() > 0)
        throw(std::runtime_error("ClusterTree::step() - currentNode_ children not empty"));
    
    vector<vector<double> > subset;
    const vector<int> &indices = currentNode_->indices;
    for(int i=0; i<indices.size(); i++) {
        subset.push_back(dataset_[indices[i]]);
    }
    currentCluster_ = new Cluster(subset, period_);
}

void ClusterTree::partitionCurrentCluster(int d) {
    currentCluster_->partition(d);
}

void ClusterTree::divideCurrentCluster() {
    vector<int> assignment = currentCluster_->cluster();
    int maxAssignmentId = *(max_element(assignment.begin(), assignment.end()));

    const vector<int> &indices = currentNode_->indices;

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
            currentNode_->children.push_back(newNode);
            queue_.push(newNode);
        }
    }

    cout << "QUEUE SIZE" << endl;
    cout << queue_.size() << endl;


    delete currentCluster_;
    currentCluster_ = NULL;
}

void ClusterTree::step() {

    if(finished()) 
        return;

    setCurrentCluster();

    if(currentCluster_->getNumPoints() < 1000) {
        delete currentCluster_;
        return;
    }

    for(int d=0; d < getNumDimensions(); d++)
        partitionCurrentCluster(d);

    // inspect partitions via getPartition and getPartitioner
    // getCurrentCluster()->getPartitioner(d) etc.

    divideCurrentCluster();
}

ClusterTree::Node& ClusterTree::getRoot() {
    return *root_;
}

Cluster& ClusterTree::getCurrentCluster() {
    return *currentCluster_;
}

int ClusterTree::getNumDimensions() const {
    return dataset_[0].size();
}