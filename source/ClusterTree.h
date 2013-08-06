#include <vector>

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

class Node {
    int population;
    std::vector<int> points;
    std::vector<Node*> children;
};

class ClusterTree {

public:

    ClusterTree(const std::vector<std::vector<double> > &dataset, const std::vector<double> &period);
    ~ClusterTree();

    // bfs search of the tree
    void stepBFS();

    // dfs search of the tree
    void stepDFS();

    void printTree();

private:

    std::vector<std::vector<double> > dataset_;
    std::vector<double> period_;
    // specifies the cluster each point belongs to
    std::vector<int> assignment_;
    // rootNode in the BFS tree
    Node *rootNode_;

    std::vector<bool> skip_;

    int counter_;
};

}