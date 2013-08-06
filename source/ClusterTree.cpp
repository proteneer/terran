#include "ClusterTree.h"
#include "Cluster.h"
#include <map>
#include <iostream>
#include <utility>

using namespace std;
using namespace Terran;

ClusterTree::ClusterTree(const vector<vector<double> > &dataset, const vector<double> &period) : 
    dataset_(dataset),
    period_(period),
    assignment_(dataset.size(), 0),
    skip_(dataset.size(), 0),
    counter_(0) {
    
}

ClusterTree::~ClusterTree() {

}

// go down one level
void ClusterTree::stepBFS() {
    
/*    
    for(int i=0; i < assignment_.size(); i++) {
        assignment_[i] = (rand() % 4)*3;
    }
*/    

    // maps assignment Id to index In Groups
    // index,           point                    period
    map <int, pair<vector<vector<double> >, vector<double> > > groups;
        
    for(int i=0; i < assignment_.size(); i++) {
       
        if(!skip_[i]) {
            (groups[assignment_[i]].first).push_back(dataset_[i]);
            (groups[assignment_[i]].second).push_back(period_[i]);
        }
    }

    for(map <int, pair<vector<vector<double> >, vector<double> > >::const_iterator it = groups.begin(); it != groups.end(); ++it) {
        cout << it->first << endl;
        const vector<vector<double> > &data = (it->second).first;
        const vector<double> &period = (it->second).second;
        
        Cluster cc(data, period);
        
        for(int d=0; d < cc.getNumDimensions(); d++) {
            cout << "dimension " << d << endl;
            cc.optimizeParameters(d);
            cc.partition(d, 0.05);

            vector<double> partitions = cc.getPartitions(d);    

            for(int i=0; i < partitions.size(); i++) {
                cout << partitions[i] << endl;
            }
            
            //testPartitions.push_back(partitions);
        }
        
        vector<int> stepAssign = cc.cluster();


    }
}