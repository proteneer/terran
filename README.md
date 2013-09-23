<h1>Terran</h1>

Terran is a clustering algorithm aimed at high dimensional spaces for continuous datasets. In particular, it will suited for datasets whose marginal distributions have distinct peaks and steep valleys. The clustering algorithm works on both periodic, aperiodic dimensions, and a mixture of both. Note that aperiodic dimensions run about ~10x slower than periodic dimension due to differences the underlying basis functions. 

The main output of the algorithm is a hierarchical cluster tree, whose leaves represent the clusters found so far. 

The algorithm has the following properties:

1) It will not oversplit dense clusters, though it certainly may undersplit.  
2) It runs in O(kD+N) time, where k is the instrinsic number of clusters, D is the number of dimensions.  

<h2>Methodology</h2>

Given a collection of points in D dimensional space, Terran first marginalizes each dimension, takes a constant subset of the resulting marginalized data, and fits a gaussian like mixture model (GMM) to the data via the Expectation Maximization algorithm. Once a sufficiently decent model is obtained, the minima of the model is found. Of the minima, it makes cuts depending on the value of each minimum: if the value is less than some sufficiently low cutoff, then a cut is made. As such, each dimension is cut into partitions. The D dimensional space is essentially divided into axis-aligned hypercubes. Each point is then assigned to the appropriate hypercube.  

<h2>Requirements</h2>

-CMake and a compatible C compiler.

Terran is named after Yutong's favorite Starcraft race.
