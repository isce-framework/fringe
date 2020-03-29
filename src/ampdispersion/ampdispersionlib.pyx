from libcpp.string cimport string

cdef extern from "ampdispersion.hpp":
    cdef cppclass ampdispersionOptions:
        ampdispersionOptions() except +
        void print()

        string inputDS
        string meanampDS
        string daDS

        int blocksize
        int memsize
        int refband

cdef extern from "ampdispersion.cpp":
    int ampdispersion_process(ampdispersionOptions*) nogil


cdef class Ampdispersion:
    '''
    Python wrapper for ampdispersion.
    '''

    cdef ampdispersionOptions *thisptr

    def __cinit__(self):
        self.thisptr = new ampdispersionOptions()


    def __dealloc__(self):
        del self.thisptr

    @property
    def inputDS(self):
        return self.thisptr.inputDS.decode('utf-8')

    @inputDS.setter
    def inputDS(self,x):
        self.thisptr.inputDS = x.encode('utf-8')


    @property
    def meanampDS(self):
        return self.thisptr.meanampDS.decode('utf-8')

    @meanampDS.setter
    def meanampDS(self,x):
        self.thisptr.meanampDS = x.encode('utf-8')

    @property
    def outputDS(self):
        return self.thisptr.daDS.decode('utf-8')

    @outputDS.setter
    def outputDS(self,x):
        self.thisptr.daDS = x.encode('utf-8')

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
    def refband(self):
        return self.thisptr.refband

    @refband.setter
    def refband(self,x):
        self.thisptr.refband = x 

    def print(self):
        self.thisptr.print()

    def run(self):
        with nogil:
            ampdispersion_process(self.thisptr)
