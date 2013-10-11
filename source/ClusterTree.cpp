#include "ClusterTree.h"
#include "Cluster.h"
#include <iostream>
#include <utility>
#include <algorithm>
#include <fstream>
#include <stdexcept>

using namespace std;
using namespace Terran;

ClusterTree::ClusterTree(const vector<vector<double> > &dataset, const vector<int> &period) : 
    dataset_(dataset),
    period_(period),
    root_(NULL),
	minPointsPerCluster_(3000),
    currentCluster_(NULL) {

	for(int i=0; i < period.size(); i++) {
		if(period[i] != 1 && period[i] != 0) {
			
			throw(std::runtime_error("period must either be zero (false), or one (true)"));
		}
	}

    vector<int> points(dataset.size());
    for(int i=0; i<points.size(); i++) {
        points[i]=i;
    }
    root_ = new Node;
    root_->indices = points;
    queue_.push(root_);

	setCurrentCluster();
}

ClusterTree::~ClusterTree() {
    delete root_;
    delete currentCluster_;
}

int ClusterTree::getNumPoints() const {
    return dataset_.size();
}

// traverse down the to the leaves
vector<int> ClusterTree::getAssignment() const {

    vector<vector<int> > clusters;
    vector<const Node*> leaves = getLeaves();

    for(int i=0; i < leaves.size(); i++) {
        clusters.push_back(leaves[i]->indices);
    }
    // check that all points exist
    vector<int> allPoints;
    for(int i=0; i < clusters.size(); i++) {
        allPoints.insert(allPoints.end(),clusters[i].begin(),  clusters[i].end());
    }
    if(allPoints.size() != getNumPoints()) {
        throw(std::runtime_error("Wrong number of points!"));
    }
    sort(allPoints.begin(), allPoints.end());

    for(int i=1; i < allPoints.size(); i++) {
        if(allPoints[i] != allPoints[i-1]+1) {
            throw(std::runtime_error("Bad point found!"));
        }
    }

    vector<int> assignment(getNumPoints());
    for(int i=0; i<clusters.size(); i++) {
        for(int j=0; j < clusters[i].size(); j++) {
            assignment[clusters[i][j]] = i;
        }
    }
    
    return assignment;
}

vector<double> ClusterTree::getPoint(int n) const {
    return dataset_[n];
}

int ClusterTree::getNumClusters() const {
    return getLeaves().size();
}

void ClusterTree::setCurrentCluster() {

    if(queue_.size() > 0) {
		if(currentCluster_ != NULL)
			throw(std::runtime_error("ClusterTree::currentCluster_ is not set to NULL, has divideCluster() been called?"));
		currentNode_ = queue_.front();
    
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
}

int ClusterTree::queueSize() const {
	return queue_.size();
}

/*
void ClusterTree::partitionCurrentCluster(int d) {
    currentCluster_->partition(d);
}
*/

void ClusterTree::divideCurrentCluster() {
	// get assignment of points from current cluster
    vector<int> assignment = currentCluster_->assign();
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
            currentNode_->children.push_back(newNode);
			
			// add this new cluster to the todo list if it has more than 3500 points.
			// (if there are less than 3500 points its hard to marginalize)
			if(subsetIndices[j].size() > minPointsPerCluster_) {
				queue_.push(newNode);
			}
        }
    }
    delete currentCluster_;
    currentCluster_ = NULL;
	// pop the queue
	queue_.pop();
	// set next cluster
	setCurrentCluster();
}

std::vector<const ClusterTree::Node*> ClusterTree::getLeaves() const {
    queue<const Node*> bfsQueue;
    bfsQueue.push(root_);
    vector<const Node*> leaves;
    while(bfsQueue.size() != 0) {
        const Node* current = bfsQueue.front();
        bfsQueue.pop();

        if(current->children.size() == 0) {
            leaves.push_back(current);
        } else {
            for(int i=0; i<current->children.size(); i++) {
                const Node* child = current->children[i];
                bfsQueue.push(child);
            }
        }
    }
    return leaves;
}

void ClusterTree::step() {
    if(queue_.size() == 0) 
        return;
	for(int d=0; d < getNumDimensions(); d++) {
		currentCluster_->partition(d);
	}
	divideCurrentCluster();
}


ClusterTree::Node& ClusterTree::getRoot() {
    return *root_;
}

Cluster& ClusterTree::getCurrentCluster() {
	if(currentCluster_ == NULL) {
		throw(std::runtime_error("ClusterTree::getCurrentCluster() - currentCluster_ is NULL"));
	}
    return *currentCluster_;
}

int ClusterTree::getNumDimensions() const {
    return dataset_[0].size();
}