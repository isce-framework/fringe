from libcpp.string cimport string

cdef extern from "phase_link.hpp":
    cdef cppclass evdOptions:
        evdOptions() except +
        void print()

        string inputDS
        string wtsDS
        string outputFolder
        string outputCompressedSlcFolder
        string compSlc
        string coherence

        int blocksize
        int memsize

        double prob

        int Nx
        int Ny
        int minNeighbors

        string method 

        int miniStackCount
        int bandWidth

cdef extern from "phase_link.cpp":
    int evd_process(evdOptions*) nogil


cdef class Phaselink:
    '''
    Python wrapper for evd.
    '''

    cdef evdOptions *thisptr

    def __cinit__(self):
        self.thisptr = new evdOptions()


    def __dealloc__(self):
        del self.thisptr

    @property
    def inputDS(self):
        return self.thisptr.inputDS.decode('utf-8')

    @inputDS.setter
    def inputDS(self,x):
        self.thisptr.inputDS = x.encode('utf-8')


    @property
    def outputFolder(self):
        return self.thisptr.outputFolder.decode('utf-8')

    @outputFolder.setter
    def outputFolder(self,x):
        self.thisptr.outputFolder = x.encode('utf-8')

    @property
    def outputCompressedSlcFolder(self):
        return self.thisptr.outputCompressedSlcFolder.decode('utf-8')

    @outputCompressedSlcFolder.setter
    def outputCompressedSlcFolder(self,x):
        self.thisptr.outputCompressedSlcFolder = x.encode('utf-8')
    
    @property
    def compSlc(self):
        return self.thisptr.compSlc.decode('utf-8')

    @compSlc.setter
    def compSlc(self,x):
        self.thisptr.compSlc = x.encode('utf-8')

    @property
    def weightsDS(self):
        return self.thisptr.wtsDS.decode('utf-8')

    @weightsDS.setter
    def weightsDS(self,x):
        self.thisptr.wtsDS = x.encode('utf-8')

    @property
    def minimumNeighbors(self):
        return self.thisptr.minNeighbors

    @minimumNeighbors.setter
    def minimumNeighbors(self,x):
        self.thisptr.minNeighbors = x

    @property
    def miniStackCount(self):
        return self.thisptr.miniStackCount

    @miniStackCount.setter
    def miniStackCount(self,x):
        self.thisptr.miniStackCount = x

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
    def method(self):
        return self.thisptr.method.decode('utf-8')

    @method.setter
    def method(self, x):
        self.thisptr.method = x.encode('utf-8')

    @property
    def bandWidth(self):
        return self.thisptr.bandWidth

    @bandWidth.setter
    def bandWidth(self,x):
        self.thisptr.bandWidth = x

    def print(self):
        self.thisptr.print()

    def run(self):
        with nogil:
            evd_process(self.thisptr)
