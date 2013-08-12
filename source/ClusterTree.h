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
    // currently found in the tree

    int getNumClusters() const;

    // return point n
    std::vector<double> getPoint(int n) const;

    bool finished() const;

    // get the set of points step() will operate on
    /*
    std::vector<std::vector<double> > getStepPoints() const;
    */

    // get the current assignment of points into clusters
    std::vector<int> getAssignment() const;

    Node& getRoot();

    Cluster& getCurrentCluster();

    // take one step in BFS
    // returns true if successful
    // returns false otherwise
    // can be easily used by a GUI system
    // 1. Given a set of points, find Partitions
    // 2. Given partitions, generate new sets of points

    // for complete automation if no fine tuning is desired.
    void step();

    // pattern:
    // call partition to inspect the dividers for each dimensions
    // check the partitions for each dimension to see how well behaved they are
    // divide the cluster into new smaller clusters
    
    // pops the queue and assigns current cluster, this removes the head element
    void setCurrentCluster();

    // partition the dth dimension in the cluster
    void partitionCurrentCluster(int d);
    
    void divideCurrentCluster();

    
private:

    // N x D
    std::vector<std::vector<double> > dataset_; 
    std::vector<double> period_;
    std::queue<Node*> queue_;
    Node* root_;

    Cluster* currentCluster_;
    Node* currentNode_;
    
};

}
