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
        
cdef class PyCluster:
    
    cdef Cluster *thisptr
    
    def __cinit__(self, vector[vector[double]] data, vector[int] period):
        self.thisptr = new Cluster(data,period)

    def __dealloc__(self):
        del self.thisptr
     
    def partition(self, int d):
        self.thisptr.partition(d)
        return self.thisptr.getPartition(d)
    
    def assign(self):
        return self.thisptr.assign()

    @property
    def dimensions(self):
        d = 0
        while d < self.thisptr.getNumDimensions():
            yield d
            d += 1

    @property  
    def shape(self):
        return self.thisptr.getNumPoints(), self.thisptr.getNumDimensions()

#cdef extern from "../include/EMGaussian.h" namespace "Terran":        
#        
#cdef extern from "../include/EM.h" namespace "Terran": 
#    cdef cppclass EM:
#        EM(vector[vector[double]]) except +
#        void setParameters(vector[Param])
#        vector[Param] getParameters()
#        void setMaxSteps(int maxSteps