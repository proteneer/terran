<h1>Terran</h1>

Terran is a clustering algorithm aimed at high dimensional spaces for continuous datasets. In particular, it will suited for datasets whose marginal distributions have distinct peaks and steep valleys. The clustering algorithm works on both periodic, aperiodic dimensions, and a mixture of both. Note that aperiodic dimensions run about ~10x slower than periodic dimension due to differences the underlying basis functions. 

The main output of the algorithm is a hierarchical cluster tree, whose leaves represent the clusters found so far. 

The algorithm has the following properties:

1) It will not oversplit dense clusters, though it certainly may undersplit.  
2) It runs in O(kD+N) time, where k is the instrinsic number of clusters, D is the number of dimensions.  

<h2>Requirements</h2>

-CMake and a compatible C compiler.

Terran is named after Yutong's favorite Starcraft race.
