# distutils: language = c++
# on windows, execute "set VS90COMNTOOLS=%VS100COMNTOOLS%"

# classes and methods exposed:
# EM
# EMGaussian
# EMPeriodicGaussian
# Cluster

from libcpp cimport bool
from libcpp.vector cimport vector
        
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

    # these are instance variables
    cdef Cluster* __thisptr
    cdef int __delete
    
    def __cinit__(self, *args):
        """
        Initializes a Cluster object.

        Arguments:
        Arg[0] -- a coordinates array of shape (N,D) where,
                  N is the number of points, D is the number dimensions
        Arg[1] -- a 1D boolean array of size D where,
                  1 denotes the domain is periodic, and 0 denotes not periodic
        """
        if len(args) == 2:
            data = <vector[vector[double]]?> args[0]
            period = <vector[int]?> args[1]
            self.__thisptr = new Cluster(data, period)
            self.__delete = 1
        
    cdef __initFromRawPointer(self, Cluster* cptr):
        self.__thisptr = cptr
        self.__delete = 0

    def __dealloc__(self):
        if self.__delete:
            del self.__thisptr
            
    def partition(self, int d):
        """
        Find deep valleys in marginalized distribution of dimension d.

        Returns a list of cutting point(s)
        """
        
        self.__thisptr.partition(d)
        return self.__thisptr.getPartition(d)
    
    def assign(self):
        """
        Assign points to a given cluster (0-indexed). 

        Must be called after partition(d) has been invoked on each dimension
        """
        return self.__thisptr.assign()
                
    property subsample:
        """ 
        Mutable: the number of points to subsample when partitioning each dimension
        """
        def __get__(self): return self.__thisptr.getSubsampleCount()
        def __set__(self, int s): self.__thisptr.setSubsampleCount(s)
    
    property shape:
        """
        Immutable: returns the shape of the underlying data
        """
        def __get__(self): return self.__thisptr.getNumPoints(), self.__thisptr.getNumDimensions()
    
    @property
    def dimensions(self):
        """
        A generator used to loop over the dimensions
        """
        d = 0
        while d < self.__thisptr.getNumDimensions():
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
    
    cdef ClusterTree *__thisptr
    
    def __cinit__(self, vector[vector[double]] data, vector[int] period):
        """
        Initializes a Cluster object.

        Arguments:
        data   -- a coordinates array of shape (N,D) where,
                  N is the number of points, D is the number dimensions
        period -- a 1D boolean array of size D where,
                  1 denotes the domain is periodic, and 0 denotes not periodic
        """
        self.__thisptr = new ClusterTree(data,period)
                
    def __dealloc__(self):
        del self.__thisptr
              
    def get_current_cluster(self):
        """
        Returns a reference to the cluster object currently being investigated.

        Note: The ownership of the object returned by this function belongs to the ClusterTree class.
              You cannot take ownership of it and/or try to delete it.
        """
        pyc = PyCluster()
        cdef Cluster* ptr = &(self.__thisptr.getCurrentCluster())
        pyc.__initFromRawPointer(ptr)
        return pyc
        