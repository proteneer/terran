#include "ClusterTree.h"
#include "Cluster.h"
#include <iostream>
#include <utility>
#include <algorithm>
#include <fstream>
#include <stdexcept>
#include <complex>

using namespace std;
using namespace Terran;

ClusterTree::ClusterTree(const vector<vector<double> > &dataset, const vector<int> &period) : 
    dataset_(dataset),
    period_(period),
    root_(NULL),
    currentCluster_(NULL),
	lastCalledFunction_(NONE) {

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
}

ClusterTree::~ClusterTree() {
    delete root_;
	delete currentCluster_;
}

int ClusterTree::getNumPoints() const {
    return dataset_.size();
}


static bool vvecSortByDescendingSize(const vector<int> &v1, const vector<int> &v2) {
	return v1.size() > v2.size();
}

vector<vector<double> > ClusterTree::centroid(int k) {

	vector<int> assignment = assign();
	int maxAssignmentId = *(max_element(assignment.begin(), assignment.end()));
	int numClusters = maxAssignmentId+1;

	if(k > numClusters) {
		throw(std::runtime_error("ClusterTree::centroid() - k > numClusters"));
	}

	vector<vector<int> > groups(numClusters);

	for(int i=0; i < dataset_.size(); i++) {
		groups[assignment[i]].push_back(i);
	}
	sort(groups.begin(), groups.end(), vvecSortByDescendingSize);
	groups.resize(k);
	
	vector<vector<double> > centroids(groups.size());

	for(int i=0; i < groups.size(); i++) {
		vector<double> center(period_.size());
		const vector<int> &points = groups[i];
		// for each dimension
		for(int d = 0; d < period_.size(); d++) {
			double mean = 0;
			// if the dimension is periodic we use directional statistics
			// to compute the periodic mean
			if(period_[d]) {
				double real = 0;
				double imag = 0;
				for(int n = 0; n < points.size(); n++) {
					double coord = dataset_[points[n]][d];
					real += cos(coord);
					imag += sin(coord);
				}
				complex<double> z(real/points.size(), imag/points.size());
				mean = arg(z);
			} else {
				for(int n = 0; n < points.size(); n++) {
					double coord = dataset_[points[n]][d];
					mean += coord;
				}
				mean = mean / points.size();
			}
			center[d] = mean;
		}
		centroids[i] = center;
	}

	return centroids;
}

// traverse down the to the leaves
vector<int> ClusterTree::assign() const {

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

int ClusterTree::getNumClusters() const {
    return getLeaves().size();
}

// setCurrentCluster must be called exactly once. 

void ClusterTree::setCurrentCluster(Partitioner* partitioner) {

	if(lastCalledFunction_ == SET_CURRENT_CLUSTER) {
		throw(std::runtime_error("ClusterTree::setCurrentCluster - called lastCalledFunction twice"));
	}

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

        if(partitioner != NULL) 
            currentCluster_ = new Cluster(subset, period_, partitioner);
        else
            currentCluster_ = new Cluster(subset, period_);
        
    }

	lastCalledFunction_ = SET_CURRENT_CLUSTER;

}

int ClusterTree::queueSize() const {
	return queue_.size();
}

void ClusterTree::divideCurrentCluster(int count) {

	if(lastCalledFunction_ == DIVIDE_CURRENT_CLUSTER) {
		throw(std::runtime_error("ClusterTree::setCurrentCluster - called divideCurrenCluster twice"));
	}

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
			if(subsetIndices[j].size() > count) {
				queue_.push(newNode);
			}
        }
    }
    delete currentCluster_;
    currentCluster_ = NULL;
	// pop the queue
	queue_.pop();

	lastCalledFunction_ = DIVIDE_CURRENT_CLUSTER;

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
    setCurrentCluster();
    if(queue_.size() == 0) 
        return;
	for(int d=0; d < getNumDimensions(); d++) {
		currentCluster_->partition(d);
	}
	divideCurrentCluster(3000);
}


ClusterTree::Node& ClusterTree::getRoot() {
    return *root_;
}

Cluster& ClusterTree::getCurrentCluster() {
	if(lastCalledFunction_ != SET_CURRENT_CLUSTER) {
		throw(std::runtime_error("ClusterTree::getCurrentCluster() - setCurrentCluster() not called yet!"));
	}
	if(currentCluster_ == NULL) {
		throw(std::runtime_error("ClusterTree::getCurrentCluster() - currentCluster_ is NULL"));
	}
    return *currentCluster_;
}

int ClusterTree::getNumDimensions() const {
    return dataset_[0].size();
}