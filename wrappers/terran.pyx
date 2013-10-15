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
    cdef int __autoDelete
    
    def __init__(self):
        """
        Do not invoke publicly
        """
        pass
       
    def __initFromArray(self, vector[vector[double]] data, vector[int] period):
        self.__thisptr = new Cluster(data, period)
        self.__autoDelete = 1
        
    cdef __initFromRawPointer(self, Cluster *ptr):
        self.__thisptr = ptr
        self.__autoDelete = 0
        
    @classmethod
    def fromArray(cls, vector[vector[double]] data, vector[int] period):
        """
        data - shape (N,D) where N is the number of points, D is the number dimensions
        period - array of size D, allowed values of period[d] are 0 and 1, where:
                                  if period[d] is 0 then dimension d is not periodic
                                  if period[d] is 1 then dimension d is periodic 
        """
        instance = cls()
        instance.__initFromArray(data, period)
        return instance
  
    def __dealloc__(self):
        if self.__autoDelete:
            del self.__thisptr
            
    def partition(self, int d):
        """
        Find deep valleys in marginalized distribution of dimension d and
        returns the cutting point(s)
        """
        
        self.__thisptr.partition(d)
        return self.__thisptr.getPartition(d)
    
    def assign(self):
        """
        Assign points to a given cluster (0-indexed). 
        Must be called after partition(d) has been invoked on each dimension
        """
        return self.__thisptr.assign()
           
    def set_subsample(self, int count):
        """
        Set the number of points to subsample when partitioning each dimension
        """
        self.__thisptr.setSubsampleCount(count)
        
    property subsample:
        def __get__(self): return self.__thisptr.getSubsampleCount()
        def __set__(self, int s): self.__thisptr.setSubsampleCount(s)
    
    property shape:
        def __get__(self): return self.__thisptr.getNumPoints(), self.__thisptr.getNumDimensions()
    
    @property
    def dimensions(self):
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
        data - shape (N,D) where N is the number of points, D is the number dimensions
        period - array of size D, allowed values of period[d] are 0 and 1, where:
                                  if period[d] is 0 then dimension d is not periodic
                                  if period[d] is 1 then dimension d is periodic 
        """
        self.__thisptr = new ClusterTree(data,period)
                
    def __dealloc__(self):
        del self.__thisptr
        
    def __gcc(self):
        pyc = PyCluster()
        cdef Cluster* ptr = &(self.__thisptr.getCurrentCluster())
        pyc.__initFromRawPointer(ptr)
        return pyc
        
    def get_current_cluster(self):
        return self.__gcc()
        