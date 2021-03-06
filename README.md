<h1>Terran</h1>

Terran is a clustering algorithm on high dimensional continuous datasets. It is well suited for datasets whose marginal distributions have distinct peaks and valleys. Each dimensions in the dataset may be either periodic or aperiodic. Note that aperiodic dimensions run about ~10x slower than periodic dimensions due to differences in the underlying basis functions.

The output of the algorithm is a hierarchical cluster tree, whose leaves represent the clusters found so far.

The algorithm has the following properties:

1) It will not oversplit dense clusters, though it certainly may undersplit.  
2) It runs in O(k(D+N)) time, where k is the instrinsic number of clusters, D is the number of dimensions.

<h2>Requirements</h2>

Terran has no dependencies on external libraries. CMake is needed to prepare the respective build environments for Windows, Linux, and OS X. Installion is via:

``` bash
$> mkdir build; cd build;
$> cmake -D CMAKE_BUILD_TYPE=RELEASE ../
$> make -j4
$> make test
```
Python wrappers can also be built via cython and distutils. 
``` bash
$> cd wrappers
$> python setup.py build_ext
```

<h2>Methodology</h2>

Given a collection of points in D dimensional space, Terran first marginalizes each dimension, takes a constant subset (about 2500 points) of the resulting marginalized data, and fits a gaussian like mixture model (GMM) to the data via the Expectation Maximization algorithm. Once a sufficient model is obtained, the minima of the model is found. Of the minima, it makes cuts depending on the value of each minimum.

As a result, each dimension is partitioned into disjoint intervals and the resulting D dimensional space is divided into axis-aligned hypercubes. Each point is then assigned to the appropriate hypercube, which then becomes its own cluster. The whole process repeats until each dimension in each cluster can no longer be partitioned.

<h3> Visual Example </h3>

Consider the following 2D dataset with 3 clusters. The marginalized distribution for each dimension is shown at the top and right of plot. An EM run on the marginalized data for dimension 0 results in a cut at _x_. Dimension 1 has no suitable cutting points.

```     
        /-\          _                                     /-\          _         
  _____/   \________/ \______                        _____/   \___x____/ \______  
 PI--------------------------PI                     PI--------------------------PI 
 |       *                   |\                     |       *     |             |\        
 |      ***                  | \                    |      ***    |             | \   
 |     **0**                 |  \                   |     **0**   c             |  \ 
 |      ***         *        |  |                   |      ***    u    *        |  | 
d|       *         ***       |   \      run EM     d|       *     t   ***       |   \
i|                **2**      |   |  -------------> i|             |  **2**      |   |
m|       *         ***       |   /                 m|       *     h   ***       |   /
 |      ***         *        |  |                   |      ***    e    *        |  |
1|     **1**                 |  /                  1|     **1**   r             |  /
 |      ***                  | /                    |      ***    e             | /
 |       *                   |/                     |       *     |             |/
-PI--------------------------PI                    -PI--------------------------PI  
             dim 0                                               dim 0  
```

The space is partitioned into two parts, note that the left plot can be further cut. The right plot is finished as no dimension can be disjointly partitioned

```     
        /-\                                                            /\         
  _____/   \_________________                        _________________/  \______  
 PI--------------------------PI                     PI--------------------------PI 
 |       *                   |\                     |                           ||
 |      ***                  | \                    |  finished                 ||  
 |     **0**                 |  |                   |                           ||  
 |      ***                  | /                    |                  *        | \
d|       *                   |/                    d|                 ***       |  \
i|---------cut here----------|x          +         i|                **2**      |   | 
m|       *                   |\                    m|                 ***       |  / 
 |      ***                  | \                    |                  *        | / 
1|     **1**                 |  |                  1|                           || 
 |      ***                  | /                    |                           || 
 |       *                   |/                     |                           ||
-PI--------------------------PI                    -PI--------------------------PI
         dim 0                                                   dim 0  
```

The left dataset is subsequently partitioned again:

```     
        /-\                                                /-\                   
  _____/   \_________________                        _____/   \_________________  
 PI--------------------------PI                     PI--------------------------PI 
 |       *                   |\                     |                           | 
 |      ***      finished    | \                    |     finished              | 
 |     **0**                 | |                    |                           |
 |      ***                  | /                    |                           | 
d|       *                   |/                    d|                           |
i|                           |                     i|                           |
m|                           |                     m|       *                   |\
 |                           |                      |      ***                  | \
1|                           |                     1|     **1**                 |  |
 |                           |                      |      ***                  | / 
 |                           |                      |       *                   |/
-PI--------------------------PI                    -PI--------------------------PI
             dim 0                                               dim 0  
```

The output is a hierarchical tree, whose leaves represent the found clusters.

<h3> C++ Example </h3>
```cpp
#include <Terran.h>

using namespace Terran;

void example() {

    // initialize input data yourself, as data is in NxD format
    vector<vector<double> > dataset;
    
    // set both dimensions to be periodic
    vector<double> periodset;
    periodset.push_back(2*PI);
    periodset.push_back(2*PI);
    
    // construct a clusterTree
    ClusterTree ct(dataset, periodset);

    // each iteration pops the BFS queue, partitions each dimension,
    // and divides the data.
    while(!ct.finished()) {
        ct.setCurrentCluster();
        ct.getCurrentCluster().partitionAll();
        ct.divideCurrentCluster(2000);
    }

    // find the point assignment
    vector<int> assignment = ct.getAssignment();
}
```
<h3> Python Example </h3>
``` python
import sys
import terran
import numpy as np
import random

f = np.load('/home/yutong/shaw_ww_dihedrals.npz', 'r')

phi = f['phis']
psi = f['psis']

phis = np.loadtxt('phis.txt', dtype=float)
psis = np.loadtxt('psis.txt', dtype=float)

dihedrals = np.concatenate((phis,psis), axis=1)

periods = np.ones(dihedrals.shape[1],dtype='int')
tree = terran.PyClusterTree(dihedrals, periods)

count = 0
while(tree.queue_size > 0):
    print("queue size: ", tree.queue_size)
    print("clusters found: ", max(tree.assign())+1)
    cluster = tree.get_current_cluster()
    cluster.partition_all()
    tree.divide_current_cluster(2000)
    assignment = tree.assign()
    # write best assignment so far
    np.savetxt('assign_loop_'+str(count)+'.txt', assignment, fmt='%d')
``` 

<h2> Misc </h2>
Terran is named after Yutong's Starcraft race.
