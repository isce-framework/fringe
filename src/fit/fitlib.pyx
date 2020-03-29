'''
Copyright (c) 2018-, California Institute of Technology ("Caltech").
U.S. Government sponsorship acknowledged.
All rights reserved.
Author (s): Heresh Fattahi
'''

import cython
import numpy as np
cimport numpy as np
cimport fitlib

np.import_array()

cdef extern from "numpy/arrayobject.h":
    void PyArray_ENABLEFLAGS (np.ndarray array, int flags)

cdef class PyFit:
    '''
    Python wrapper for Fit class
    '''
    cdef Fit* c_fit

    def __cinit__(self):
        self.c_fit = new Fit()

    def __dealloc__ (self):
        del self.c_fit

    @property
    def inputDS(self):
        return self.c_fit.inputDS.decode('utf-8')

    @inputDS.setter
    def inputDS(self,x):
        self.c_fit.inputDS = x.encode('utf-8')

    @property
    def outputDS(self):
        return self.c_fit.outputDS.decode('utf-8')

    @outputDS.setter
    def outputDS(self,x):
        self.c_fit.outputDS = x.encode('utf-8') 

    @property
    def blocksize(self):
        return self.c_fit.blocksize

    @blocksize.setter
    def blocksize(self,x):
        self.c_fit.blocksize = x 

    @property
    def earthquake_time(self):
        return self.c_fit.t_Eq

    @earthquake_time.setter
    def earthquake_time(self,x):
        self.c_fit.t_Eq = x

    @property
    def Tau(self):
        return self.c_fit.Tau

    @Tau.setter
    def Tau(self,x):
        self.c_fit.Tau = x

    def set_blocksize(self, x):
        self.c_fit.set_blocksize(x)

    def set_outputDataset(self, x):
        self.c_fit.set_outputDataset(x.encode('utf-8'))

    def set_earthquake_time(self, x):
        self.c_fit.set_earthquake_time(x)

    def set_tau(self, x):
        self.c_fit.set_Tau(x)

    def set_Npar(self, x):
        self.c_fit.set_Npar(x)

    def open_inputDataset(self, x):
        self.c_fit.OpenDataset(x.encode('utf-8'))

    def fit(self, modelType):
        self.c_fit.get_time()
        self.c_fit.CreateDataset()
        self.c_fit.fit_timeseries(modelType.encode('utf-8'))
        self.c_fit.CloseDataset()
    
    def set_G(self, np.ndarray[double, ndim=2, mode="c"] G not None):
        self.c_fit.set_G(&G[0,0], G.shape[0], G.shape[1])

    def fit_timeseries(self):
        self.c_fit.get_time()
        self.c_fit.CreateDataset()
        self.c_fit.fit_timeseries()
        self.c_fit.CloseDataset()

