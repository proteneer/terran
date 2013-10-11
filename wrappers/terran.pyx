# distutils: language = c++
# on windows, execute "set VS90COMNTOOLS=%VS100COMNTOOLS%"

# classes and methods exposed:
# EM
# EMGaussian
# EMPeriodicGaussian
# Cluster

from libcpp cimport bool
from libcpp.vector cimport vector

cdef extern from "../include/Param.h" namespace "Terran":
    cdef cppclass Param:
        Param(double, double, double) except +
        double p,u,s

# In python, the most intuitive thing is probably to just pass around 3-tuples. 
cdef class PyParam:
    cdef Param *thisptr
    def __cinit__(self, double p, double u, double s):
        self.thisptr = new Param(p,u,s)
    def __dealloc__(self):
        del self.thisptr
    
    property p:
        def __get__(self): return self.thisptr.p
        def __set__(self, p): self.thisptr.p = p
        
    property u:
        def __get__(self): return self.thisptr.u
        def __set__(self, u): self.thisptr.u = u
        
    property s:
        def __get__(self): return self.thisptr.s
        def __set__(self, s): self.thisptr.s = s
        
cdef extern from "../include/Cluster.h" namespace "Terran":
    cdef cppclass Cluster:
        Cluster(vector[vector[double]],vector[int]) except +
        int getNumDimensions()
        int getNumPoints()
        bool isPeriodic(int)
        void partition(int)
        vector[double] getPartition(int)
        vector[int] assign()
        void setSubsampleCount(int)
        int getSubsampleCount()
           
cdef class PyCluster:
    
    cdef Cluster* thisptr
        
    def __cinit__(self, vector[vector[double]] data, vector[int] period):
        """
        data - shape (N,D) where N is the number of points, D is the number dimensions
        period - array of size D, allowed values of period[d] are 0 and 1, where:
                                  if period[d] is 0 then dimension d is not periodic
                                  if period[d] is 1 then dimension d is periodic 
        """
        self.thisptr = new Cluster(data,period)
        
    def __cinit__(self):
        """
        Do not use this method. This constructor is mainly used to help out PyClusterTree
        """
        pass
        
    def __dealloc__(self):
        del self.thisptr
     
    def partition(self, int d):
        """
        Find deep valleys in marginalized distribution of dimension d and
        returns the cutting point(s)
        """
        self.thisptr.partition(d)
        return self.thisptr.getPartition(d)
    
    def assign(self):
        """
        Assign points to a given cluster (0-indexed). 
        Must be called after partition(d) has been invoked on each dimension
        """
        return self.thisptr.assign()

    def set_subsample(self, int count):
        """
        Set the number of points to subsample when partitioning each dimension
        """
        self.thisptr.setSubsampleCount(count)

    property subsample:
        def __get__(self): return self.thisptr.getSubsampleCount()
        def __set__(self, int s): self.thisptr.setSubsampleCount(s)
    
    property shape:
        def __get__(self): return self.thisptr.getNumPoints(), self.thisptr.getNumDimensions()
    
    @property
    def dimensions(self):
        d = 0
        while d < self.thisptr.getNumDimensions():
            yield d
            d += 1

cdef extern from "../include/ClusterTree.h" namespace "Terran":
    cdef cppclass ClusterTree:
        ClusterTree(vector[vector[double]],vector[int]) except +
        int getNumDimensions()
        int getNumPoints()
        int getNumClusters()
        vector[int] assign()
        Cluster& getCurrentCluster()
        
cdef class PyClusterTree:
    
    cdef ClusterTree *thisptr
    
    def __cinit__(self, vector[vector[double]] data, vector[int] period):
        """
        data - shape (N,D) where N is the number of points, D is the number dimensions
        period - array of size D, allowed values of period[d] are 0 and 1, where:
                                  if period[d] is 0 then dimension d is not periodic
                                  if period[d] is 1 then dimension d is periodic 
        """
        self.thisptr = new ClusterTree(data,period)
                
    def __dealloc__(self):
        del self.thisptr
        
    def get_current_cluster(self):
        cluster = PyCluster()
        cluster.thisptr = &self.thisptr.getCurrentCluster();
        return cluster
