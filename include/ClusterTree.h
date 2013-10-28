#include <vector>
#include <queue>

#include "export.h"
#include "Cluster.h"

/* A BFS based ClusterTree class. 

 Usage:

 ClusterTree ct(dataset, period);

 The resulting tree might look like:
                     *                      
                     0                          
           __________|_________*                
           0         1         2
                           *___|___*
                           0   1   2
                          / \     / \
                         0   1   0   1

 Each level of the BFS tree uses at most O(N*d) space.

 */
namespace Terran {

class TERRAN_EXPORT ClusterTree {

public:

    struct Node {
        Node() {};
        ~Node() {
            for(int i=0;i<children.size(); i++) {
                delete children[i];
                children[i] = NULL;
            }
                
        };
        // points that belong to this cluster
        std::vector<int> indices; 
        // partition dividers
        std::vector<std::vector<double> > partitions;
        // the clusters results from this
        std::vector<Node*> children;
    };

    // dataset: NxD matrix, N is number of points, D is dimension
    // period denotes the periodicity 
    ClusterTree(const std::vector<std::vector<double> > &dataset, const std::vector<int> &isPeriodic);

    ~ClusterTree();

	// compute the centroid of the first k most populated cluster
	/*
	std::vector<std::vector<double> > centroid(int k);
	*/

    // return number of points in the dataset
    int getNumPoints() const;

    // return number of dimensions in the dataset
    int getNumDimensions() const;

    // get number of clusters found so far, equal to number of leaves 
    int getNumClusters() const;
	
    // get the current assignment of points into clusters
    std::vector<int> assign() const;

    // get a reference to the currently processed cluster
	Cluster& getCurrentCluster();

	// return number of clusters left to process
	int queueSize() const;

    // the step function is convenience function that does the following
    // 1. see if the queue is empty, if not, pop the head and set it as the currentCluster
    // 2. if the currentCluster has less than 500 points, then mark as "finished"
    // 3. partition each dimension in currentCluster using simpleRun()
    // 4. it divides the current cluster into subcluster
    // for fine tuned control, all four steps can be reconfigured as needed
    void step();

	// divide the current cluster into more clusters by looking at the assignments
	// cannot be called twice in a row
	void divideCurrentCluster(int count);

    // set the current cluster, can pass in optional partitioner 
	// cannot be called twice in a row
    void setCurrentCluster(Partitioner* partitioner = NULL);

private:
	
	enum CalledFunction { NONE, SET_CURRENT_CLUSTER, DIVIDE_CURRENT_CLUSTER };
	CalledFunction lastCalledFunction_;

    // returns a reference to the root of the BFS tree
    Node& getRoot();

    void validateData();
	
    // returns found leaves so far
    std::vector<const Node*> getLeaves() const;

    // N x D
    std::vector<std::vector<double> > dataset_; 
    std::vector<int> period_;
    std::queue<Node*> queue_;
    Node* root_;
    Cluster* currentCluster_;
    Node* currentNode_;
    
};

}
