
import numpy as np

from scipy.cluster.vq import kmeans2

from scipy.spatial.distance import pdist
from scipy.spatial.distance import squareform
import pylab as pl
import csv
import os as os

from math import floor
from math import sqrt

from sklearn import cluster, datasets, metrics
from sklearn.metrics import euclidean_distances
from sklearn.neighbors import kneighbors_graph
from sklearn.preprocessing import StandardScaler

period = [2*3.14159265, 2*3.14159265]
periodNP = np.array([2*3.14159265, 2*3.14159265])

def loadData(filename):
    for line in csv.reader(open(filename), delimiter=' '):
        if line:
            yield line

def periodicDistance(xs, ys, ps):
    ''' given two d dimensional points x and y, of which some dimensions may be periodic, compute the
    periodic difference '''
    
    assert len(xs) == len(ys)
    assert len(ys) == len(ps)

    dist = 0;
    for (x,y,p) in zip(xs, ys, ps):
        diff = x - y
        if p != 0:
            diff -= floor(diff/p + 0.5) * p
        dist += diff*diff
    return sqrt(dist)
        
def periodicDistanceNP(xs,ys):
    diff = xs - ys
    for i in xrange(diff.shape[0]):
       if(periodNP[i] != 0):
            diff[i] -= np.floor(diff[i]/periodNP[i]+0.5)*periodNP[i]
                   
    return np.sqrt(np.dot(diff, diff))
            
def testPeriodicDistance():
    a = 0.1
    b = 0.9
    p = 1.0
    
    alist = [a]
    blist = [b]
    plist = [p]
    
    print periodicDistance(alist,blist,plist)

def computePairwiseDistances(data, period):
    distances = []
    count = 0
    for x in data:
        count += 1
        print count
        rowDistances = []
        for y in data:
            rowDistances.append(periodicDistance(x,y, period))
        distances.append(rowDistances) 
    return distances
         
def drange(start, stop, step):
    r = start
    while r < stop:
        yield r
        r += step

np.random.seed(0)
         
#read input file for terran 
highDDataset = []
for data in loadData("100d.txt"):
    highDDataset.append([float(x) for x in data])

D100 = np.array(highDDataset)


highDDatasetMulti = []
for data in loadData("100dMulti.txt"):
    highDDatasetMulti.append([float(x) for x in data])

D100Multi = np.array(highDDatasetMulti)

for idx, X in enumerate([D100Multi, D100]):

    labels_true = []
    if idx == 0:
        for assign in loadData("100dMultiLabels.txt"):
            labels_true.append(int(assign[0]))
    elif idx == 1:
        for assign in loadData("100dLabels.txt"):
            labels_true.append(int(assign[0]))    

    dbscan_low  = cluster.DBSCAN(eps=1.0, metric='euclidean')
    dbscan_high = cluster.DBSCAN(eps=2.0, metric='euclidean')
    spectral2 = cluster.SpectralClustering(n_clusters=2, 
                                          eigen_solver='arpack',
                                          affinity="nearest_neighbors")  

    spectral17 = cluster.SpectralClustering(n_clusters=17, 
                                          eigen_solver='arpack',
                                          affinity="nearest_neighbors") 

    spectral20 = cluster.SpectralClustering(n_clusters=20, 
                                          eigen_solver='arpack',
                                          affinity="nearest_neighbors")  

    spectral23 = cluster.SpectralClustering(n_clusters=23, 
                                          eigen_solver='arpack',
                                          affinity="nearest_neighbors") 

    kmeans2 = cluster.KMeans(n_clusters=2)
    kmeans17 = cluster.KMeans(n_clusters=17)
    kmeans20 = cluster.KMeans(n_clusters=20)
    kmeans23 = cluster.KMeans(n_clusters=23)

    # dbscan and spectral can easily screw up
    #for algorithm in [dbscan_low, dbscan_high, spectral, two_means, twenty_means]:

    for algorithm in [dbscan_low, dbscan_high, spectral2, spectral17, spectral20, spectral23, 
                      kmeans2, kmeans17, kmeans20, kmeans23]:
        print str(algorithm).split('(')[0]
        algorithm.fit(X)
        label = algorithm.labels_.astype(np.int)
        llist = [0]*(max(label)+1)
        if(max(label) == -1):
            print "No Clusters Found"
        else:

            print 'RI SCORE: ', metrics.adjusted_rand_score(labels_true, label)  

            #for a in label:
            #    llist[a] += 1
            #clusterIndex = 0
            #for i in llist:
            #    print 'cluster', clusterIndex, 'has', i, 'points'
            #    clusterIndex += 1

    multi_terran_labels = []
    if idx == 0:
        for assign in loadData("100dMultiTerranAssign.txt"):
            multi_terran_labels.append(int(assign[0]))
    elif idx == 1:
        for assign in loadData("100dTerranAssign.txt"):
            multi_terran_labels.append(int(assign[0]))    

    print 'TERRAN RI SCORE: ', metrics.adjusted_rand_score(labels_true, multi_terran_labels)

    #print max(label)