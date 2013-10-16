# classes and methods exposed:
# EM
# EMGaussian
# EMPeriodicGaussian
# Cluster
# ClusterTree

from libcpp cimport bool
from libcpp.vector cimport vector

# Abstract Class
cdef extern from "../include/Partitioner.h" namespace "Terran":
    cdef cppclass Partitioner:
        Partitioner() except +
        
# Concrete EM-based Partitioner class        
cdef extern from "../include/PartitionerEM.h" namespace "Terran":
    cdef cppclass PartitionerEM(Partitioner):
        PartitionerEM() except +
        void setDataAndPeriod(vector[double], bool)
        vector[double] partition()
        void evaluateModel(vector[double]&, vector[double]&, int)
        void setPartitionCutoff(double val)
        double getPartitionCutoff()
        void setInitialK(int count)
        int getInitialK()

cdef class PyPartitionerEM:
    
    cdef PartitionerEM* __thisptr
    cdef int __delete
    
    def __cinit__(self):
        self.__thisptr = new PartitionerEM()
        self.__delete = 1
    
    def __dealloc__(self):
        if self.__delete:
            del self.__thisptr 
    
    def model(self, int nsamples):
        """ Generate a curve using the underlying data model """
        cdef vector[double] xvals
        cdef vector[double] yvals
        self.__thisptr.evaluateModel(xvals, yvals, nsamples)
        return (xvals, yvals)
        
    def setDataAndPeriod(self, vector[double] data, bool periodic):
        """
        Set the data and period information.
        """
        self.__thisptr.setDataAndPeriod(data, periodic)
        
    def partition(self):
        return self.__thisptr.partition()
    
    property cutoff:
        """ 
        Mutable: minima point must be be lower than cutoff in order to be cut 
        """
        def __get__(self): return self.__thisptr.getPartitionCutoff()
        def __set__(self, double c): self.__thisptr.setPartitionCutoff(c) 
    
    property initial_k:
        """ 
        Mutable: the number of clusters to use initially when optimizing
        """
        def __get__(self): return self.__thisptr.getInitialK()
        def __set__(self, int k): self.__thisptr.setInitialK(k)
    
    
cdef extern from "../include/Cluster.h" namespace "Terran":
    cdef cppclass Cluster:
        Cluster(vector[vector[double]],vector[int]) except +
        Cluster(vector[vector[double]],vector[int], Partitioner* partitioner) except +
        int getNumDimensions()
        int getNumPoints()
        void partition(int) except+
        vector[double] getPartition(int)
        vector[int] assign() except+
        void setSubsampleCount(int)
        int getSubsampleCount()
        Partitioner& getPartitioner()
           
cdef class PyCluster:

    # Note: these are instance variables
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
        Arg[2] -- Optional: a PyPartitionerEM object. Note that PyCluster takes ownership
                  of the PyPartitionerEM object. 
        """
        
        print('Calling PyCluster() __cinit__')

        if len(args) == 2:
            data = <vector[vector[double]]?> args[0]
            period = <vector[int]?> args[1]
            self.__thisptr = new Cluster(data, period)
            self.__delete = 1
        if len(args) == 3:
            data = <vector[vector[double]]?> args[0]
            period = <vector[int]?> args[1]   
            pyPartitioner = <PyPartitionerEM?>args[2]
            pyPartitioner.__delete = 0
            self.__thisptr = new Cluster(data, period, pyPartitioner.__thisptr)
            self.__delete = 1
            
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
              
    def get_partitioner(self):
        """
        Returns a reference to the cluster object currently being investigated.

        Note: The ownership of the object returned by this function belongs to the ClusterTree class.
              You cannot take ownership of it and/or try to delete it.
        """
        #equivalent of a dynamic-cast, next in try catch for more types of partition classes
        cdef PartitionerEM *pem = <PartitionerEM*?>&(self.__thisptr.getPartitioner())
        partem = PyPartitionerEM()
        partem.__thisptr = pem
        partem.__delete = 0
        return partem
              
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
        int queueSize()
        vector[int] assign()
        Cluster& getCurrentCluster()
        void divideCurrentCluster(int)
        
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
        pyc.__thisptr = ptr
        pyc.__delete = 0
        return pyc

    def divide_current_cluster(self, int cutoff):
        """
        Divide the currently assigned cluster into more clusters, appending it into the queue.
        
        Note: cutoff is recommended to be 2000 or more, as lack of points
              makes EM and subsequent marginalization difficult
        """
        self.__thisptr.divideCurrentCluster(cutoff)

    property clusters_found:
        """
        Immutable: return number of clusters found so far
        """
        def __get__(self): return self.__thisptr.getNumClusters()

    property queue_size:
        """
        Immutable: return number of clusters to be processed in the queue
        """
        def __get__(self): return self.__thisptr.queueSize()

    property shape:
        """
        Immutable: returns the shape of the underlying data
        """
        def __get__(self): return self.__thisptr.getNumPoints(), self.__thisptr.getNumDimensions()   

    def assign(self):
        """
        Assign points to a given cluster (0-indexed) by scanning the leaves of the ClusterTree 
        """
        return self.__thisptr.assign()
        