from libcpp.string cimport string
from libcpp cimport bool

cdef extern from "despeck.hpp":
    cdef cppclass despeckOptions:
        despeckOptions() except +
        void print()

        string inputDS
        string wtsDS
        string outputDS
       
        int ibands[2]
        int blocksize
        int memsize

        int Nx
        int Ny

        bool computeCoherence

cdef extern from "despeck.cpp":
    int despeck_process(despeckOptions*) nogil


cdef class Despeck:
    '''
    Python wrapper for despeck.
    '''

    cdef despeckOptions *thisptr

    def __cinit__(self):
        self.thisptr = new despeckOptions()


    def __dealloc__(self):
        del self.thisptr

    @property
    def inputDS(self):
        return self.thisptr.inputDS.decode('utf-8')

    @inputDS.setter
    def inputDS(self,x):
        self.thisptr.inputDS = x.encode('utf-8')

    @property
    def weightsDS(self):
        return self.thisptr.wtsDS.decode('utf-8')

    @weightsDS.setter
    def weightsDS(self,x):
        self.thisptr.wtsDS = x.encode('utf-8')

    @property
    def outputDS(self):
        return self.thisptr.outputDS.decode('utf-8')

    @outputDS.setter
    def outputDS(self,x):
        self.thisptr.outputDS = x.encode('utf-8')

    @property
    def blocksize(self):
        return self.thisptr.blocksize

    @blocksize.setter
    def blocksize(self,x):
        self.thisptr.blocksize = x 

    @property
    def memsize(self):
        return self.thisptr.memsize

    @memsize.setter
    def memsize(self,x):
        self.thisptr.memsize = x

    @property
    def halfWindowX(self):
        return self.thisptr.Nx

    @halfWindowX.setter
    def halfWindowX(self,x):
        self.thisptr.Nx = x

    @property
    def halfWindowY(self):
        return self.thisptr.Ny

    @halfWindowY.setter
    def halfWindowY(self,x):
        self.thisptr.Ny = x

    
    @property
    def band1(self):
        return self.thisptr.ibands[0]

    @band1.setter
    def band1(self,x):
        self.thisptr.ibands[0] = x


    @property
    def band2(self):
        return self.thisptr.ibands[1]

    @band2.setter
    def band2(self,x):
        self.thisptr.ibands[1] = x 

    @property
    def coherenceFlag(self):
        return self.thisptr.computeCoherence

    @coherenceFlag.setter
    def coherenceFlag(self, x):
        self.thisptr.computeCoherence = x > 0

    def print(self):
        self.thisptr.print()

    def run(self):
        with nogil:
            despeck_process(self.thisptr)
