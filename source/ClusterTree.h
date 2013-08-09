#include <vector>
#include <queue>

// A BFS based hierarchical clustering algorithm

// The resulting tree might look like:
//                     *                      
//                     0                          
//           __________|_________*                
//           0         1         2
//                           *___|___*
//                           0   1   2
//                          / \     / \
//                         0   1   0   1
//
// Each level of the BFS tree uses at most O(N*d) space.

namespace Terran {

class ClusterTree {

public:

    struct ClusterTreeNode {
        ~ClusterTreeNode() {
            for(int i=0;i<children.size(); i++)
                delete children[i];
        }
        // points that belong to this cluster
        std::vector<int> indices; 
        // partition dividers
        std::vector<std::vector<double> > partitions;
        // the clusters results from this
        std::vector<ClusterTreeNode*> children;
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

    // std::vector<int> getPointsInCluster(int clusterIndex) const;

    // int getPopulation(int clusterIndex);

    // get the current assignment of points into clusters
    std::vector<int> getAssignment() const;

    // take one step in BFS
    // returns true if successful
    // returns false otherwise
    bool step();

    void getTree();
    
private:

    // N x D
    std::vector<std::vector<double> > dataset_;
    std::vector<double> period_;

    // Each cluster is comprised of point Ids, and a boolean skip flag where: 
    // TODO: change to enums later
    // true  - denotes the cluster can be skipped
    // false - denotes the cluster still needs to be searched
    //std::vector<std::pair<std::vector<int>, bool> > clusters_; 
    //std::vector<std::vector<int>, bool> > clusters_; 

    // used in BFS search
    // elements in the queue
    std::queue<ClusterTreeNode*> queue_;

    // the root of the tree
    ClusterTreeNode* root_;
    
};

}
