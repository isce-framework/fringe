from libcpp cimport bool as bool_t
from libcpp.string cimport string

cdef extern from "nmap.hpp":
    cdef cppclass nmapOptions:
        nmapOptions() except +
        void print()

        string inputDS
        string ncountDS
        string wtsDS
        string maskDS
        string method
        int blocksize
        int memsize

        double prob
        
        int Nx
        int Ny

        bool_t noGPU

cdef extern from "nmap.cpp":
    int nmap_process(nmapOptions*) nogil

cdef class Nmap:
    '''
    Python wrapper for nmap.
    '''

    cdef nmapOptions *thisptr

    def __cinit__(self):
        self.thisptr = new nmapOptions()


    def __dealloc__(self):
        del self.thisptr

    @property
    def inputDS(self):
        return self.thisptr.inputDS.decode('utf-8')

    @inputDS.setter
    def inputDS(self,x):
        self.thisptr.inputDS = x.encode('utf-8')


    @property
    def countDS(self):
        return self.thisptr.ncountDS.decode('utf-8')

    @countDS.setter
    def countDS(self,x):
        self.thisptr.ncountDS = x.encode('utf-8')

    @property
    def weightsDS(self):
        return self.thisptr.wtsDS.decode('utf-8')

    @weightsDS.setter
    def weightsDS(self,x):
        self.thisptr.wtsDS = x.encode('utf-8')

    @property
    def maskDS(self):
        return self.thisptr.maskDS.decode('utf-8')

    @maskDS.setter
    def maskDS(self,x):
        self.thisptr.maskDS = x.encode('utf-8')

    @property
    def method(self):
        return self.thisptr.method.decode('utf-8')

    @method.setter
    def method(self,x):
        self.thisptr.method = x.encode('utf-8')

    @property
    def minimumProbability(self):
        return self.thisptr.prob

    @minimumProbability.setter
    def minimumProbability(self,x):
        self.thisptr.prob = x

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
    def noGPU(self):
        return self.thisptr.noGPU

    @noGPU.setter
    def noGPU(self, x):
        self.thisptr.noGPU = x

    def print(self):
        self.thisptr.print()

    def run(self):
        with nogil:
            nmap_process(self.thisptr)
