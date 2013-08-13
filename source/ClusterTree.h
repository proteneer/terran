#include <vector>
#include <queue>

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

class ClusterTree {

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

    ClusterTree(const std::vector<std::vector<double> > &dataset, const std::vector<double> &period);

    ~ClusterTree();

    // return number of points in the dataset
    int getNumPoints() const;

    // return number of dimensions in the dataset
    int getNumDimensions() const;

    // get number of clusters found so far, equal to number of leaves 
    int getNumClusters() const;

    // return point n
    std::vector<double> getPoint(int n) const;

    // returns true if queue_ is empty, false otherwise
    bool finished() const;

    // get the current assignment of points into clusters
    std::vector<int> getAssignment() const;

    // returns a reference to the root of the BFS tree
    Node& getRoot();

    // returns a reference to the currently processed cluster
    Cluster& getCurrentCluster();

    // process one cluster if no fine tuning is desired.
    // step() also ignores clusters without sufficient points
    void step();
    
    // the following three methods provide fine grained control
    
    // pops the queue and assigns current cluster, this removes the head element
    void setCurrentCluster();

    // find partitions in dimension d
    void partitionCurrentCluster(int d);
    
    // divide the current cluster possibly into more clusters if there are
    // multiple points that can be assigned
    void divideCurrentCluster();

private:

    // returns found leaves so far
    std::vector<const Node*> getLeaves() const;

    // N x D
    std::vector<std::vector<double> > dataset_; 
    std::vector<double> period_;
    std::queue<Node*> queue_;
    Node* root_;

    Cluster* currentCluster_;
    Node* currentNode_;
    
};

}
