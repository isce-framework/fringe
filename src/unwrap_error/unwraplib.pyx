'''
Copyright (c) 2018-, California Institute of Technology ("Caltech").
U.S. Government sponsorship acknowledged.
All rights reserved.
Author (s): Heresh Fattahi
'''

import cython
import numpy as np
cimport numpy as np
cimport unwraplib

np.import_array()

cdef extern from "numpy/arrayobject.h":
    void PyArray_ENABLEFLAGS (np.ndarray array, int flags)

cdef class PyUnwrap:
    '''
    Python wrapper for Unwrap class
    '''
    cdef Unwrap* c_unwrap

    def __cinit__(self):
        self.c_unwrap = new Unwrap()

    def __dealloc__ (self):
        del self.c_unwrap

    @property
    def inputDS(self):
        return self.c_unwrap.connCompDS.decode('utf-8')

    @inputDS.setter
    def inputDS(self,x):
        self.c_unwrap.connCompDS = x.encode('utf-8')

    @property
    def maskDS(self):
        return self.c_unwrap.maskDS.decode('utf-8')

    @maskDS.setter
    def maskDS(self,x):
        self.c_unwrap.maskDS = x.encode('utf-8') 

    @property
    def blocksize(self):
        return self.c_unwrap.blocksize

    @blocksize.setter
    def blocksize(self,x):
        self.c_unwrap.blocksize = x 

    def Set_connComp_dataset(self, x):
        self.c_unwrap.Set_connComp_dataset(x.encode('utf-8'))

    def ComputeCommonMask(self):
        self.c_unwrap.OpenDataset()
        self.c_unwrap.CreateMaskDataset()
        self.c_unwrap.ComputeCommonMask()
        self.c_unwrap.CloseDataset()

    #def set_G(self, np.ndarray[double, ndim=2, mode="c"] G not None):
    #    self.c_unwrap.set_G(&G[0,0], G.shape[0], G.shape[1])

