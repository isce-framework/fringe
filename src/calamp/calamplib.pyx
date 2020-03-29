from libcpp cimport bool
from libcpp.string cimport string

cdef extern from "calamp.hpp":
    cdef cppclass calampOptions:
        calampOptions() except +
        void print()

        string inputDS
        string outputDS
        string maskDS
        double defaultValue
        bool applySqrt
        int blocksize
        int memsize

cdef extern from "calamp.cpp":
    int calamp_process(calampOptions*) nogil


cdef class Calamp:
    '''
    Python wrapper to calamp.
    '''

    cdef calampOptions *thisptr

    def __cinit__(self):
        self.thisptr = new calampOptions()


    def __dealloc__(self):
        del self.thisptr

    @property
    def inputDS(self):
        return self.thisptr.inputDS.decode('utf-8')

    @inputDS.setter
    def inputDS(self,x):
        self.thisptr.inputDS = x.encode('utf-8')


    @property
    def outputDS(self):
        return self.thisptr.outputDS.decode('utf-8')

    @outputDS.setter
    def outputDS(self,x):
        self.thisptr.outputDS = x.encode('utf-8')

    @property
    def maskDS(self):
        return self.thisptr.maskDS.decode('utf-8')

    @maskDS.setter
    def maskDS(self,x):
        self.thisptr.maskDS = x.encode('utf-8')

    @property
    def defaultValue(self):
        return self.thisptr.defaultValue

    @defaultValue.setter
    def defaultValue(self,x):
        self.thisptr.defaultValue = x

    @property
    def applySqrt(self):
        return self.thisptr.applySqrt

    @applySqrt.setter
    def applySqrt(self,x):
        self.thisptr.applySqrt = x

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

    def print(self):
        self.thisptr.print()

    def run(self):
        with nogil:
            calamp_process(self.thisptr)
