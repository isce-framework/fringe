from libcpp.string cimport string
from libcpp.vector cimport vector
from libcpp cimport bool as bool_t

cdef extern from "sbas.hpp":
    cdef cppclass Scene:
        string date
        string bperpName

    cdef cppclass Pair:
        Pair() except +
        void print()

        string masterDate
        string slaveDate
        string ifgName
        string cohName
        float threshold
        float scale
        float referenceOffset
        float deltaT

    cdef cppclass sbasOptions:
        sbasOptions() except +
        void print()

        int nPairs
        vector[Pair] pairs
        vector[Scene] scenes
        vector[int] refbox
        vector[int] bbox
        string outDir
        string refDate
        int nSAR
        vector[string] dates
        bool_t estimateDEMError
        double startingRange
        double rangeSpacing
        double wavelength
        string incAngleFile

cdef extern from "sbas.cpp":
    void c_process "process" (sbasOptions*) nogil


cdef class SceneWrapper:
    '''
    Python wrapper for scene.
    '''
    cdef Scene *thisptr
    cdef bool_t owner

    def __cinit__(self):
        self.thisptr = new Scene()
        self.owner = True

    def __dealloc__(self):
        if self.owner:
            del self.thisptr

    @staticmethod
    cdef bind(Scene *x):
        new_scene = SceneWrapper()
        del new_scene.thisptr
        new_scene.thisptr = x 
        new_scene.owner = False
        return new_scene

    @property
    def date(self):
        return self.thisptr.date.decode('utf-8')

    @date.setter
    def date(self, x):
        self.thisptr.date = x.encode('utf-8')

    @property
    def baselineFile(self):
        return self.thisptr.bperpName.decode('utf-8')

    @baselineFile.setter
    def baselineFile(self, x):
        self.thisptr.bperpName = x.encode('utf-8')


cdef class PairWrapper:
    '''
    Python wrapper for Pair.
    '''

    cdef Pair *thisptr
    cdef bool_t owner

    def __cinit__(self):
        self.thisptr = new Pair()
        self.owner = True

    def __dealloc__(self):
        if self.owner:
            del self.thisptr

    @staticmethod
    cdef bind(Pair* x):
        new_pair = PairWrapper()
        del new_pair.thisptr
        new_pair.thisptr = x
        new_pair.owner = False
        return new_pair

    def print(self):
        self.thisptr.print()

    @property
    def masterDate(self):
        return self.thisptr.masterDate.decode('utf-8')

    @masterDate.setter
    def masterDate(self, x):
        self.thisptr.masterDate = x.encode('utf-8')

    @property
    def slaveDate(self):
        return self.thisptr.slaveDate.decode('utf-8')

    @slaveDate.setter
    def slaveDate(self, x):
        self.thisptr.slaveDate = x.encode('utf-8')

    @property
    def ifgName(self):
        return self.thisptr.ifgName.decode('utf-8')

    @ifgName.setter
    def ifgName(self, x):
        self.thisptr.ifgName = x.encode('utf-8')

    @property
    def cohName(self):
        return self.thisptr.cohName.decode('utf-8')

    @cohName.setter
    def cohName(self, x):
        self.thisptr.cohName = x.encode('utf-8')

    @property
    def threshold(self):
        return self.thisptr.threshold

    @threshold.setter
    def threshold(self, x):
        self.thisptr.threshold = float(x)

    @property
    def scale(self):
        return self.thisptr.scale

    @scale.setter
    def scale(self, x):
        self.thisptr.scale = float(x)

    @property
    def deltaT(self):
        return self.thisptr.deltaT

    @property
    def referenceOffset(self):
        return self.thisptr.referenceOffset


cdef class SBASWrapper:
    '''
    Python wrapper for sbas.
    '''

    cdef sbasOptions *thisptr

    def __cinit__(self):
        self.thisptr = new sbasOptions()

    def __dealloc__(self):
        del self.thisptr

    def process(self):
        c_process(self.thisptr)

    def setNumberOfPairs(self, insize):
        '''
        Set size of input vectors.
        '''
        self.thisptr.nPairs = int(insize)
        self.thisptr.pairs.resize( int(insize))

    def setNumberOfScenes(self, insize):
        '''
        Set size of input vector.
        '''
        self.thisptr.scenes.resize(int(insize))

    @property
    def numberOfPairs(self):
        return self.thisptr.nPairs

    @property
    def numberOfImages(self):
        return self.thisptr.nSAR

    def getPair(self, index):
        return  PairWrapper.bind(&(self.thisptr.pairs[index]))

    def getScene(self, index):
        return SceneWrapper.bind(&(self.thisptr.scenes[index]))

    @property
    def referenceDate(self):
        return self.thisptr.refDate

    @referenceDate.setter
    def referenceDate(self, x):
        self.thisptr.refDate = x.encode('utf-8')

    @property
    def outputDir(self):
        return self.thisptr.outDir

    @outputDir.setter
    def outputDir(self, x):
        self.thisptr.outDir = x.encode('utf-8')

    @property
    def dates(self):
        return self.thisptr.dates

    @property
    def referenceBox(self):
        return self.thisptr.refbox

    @referenceBox.setter
    def referenceBox(self, x):
        if len(x) == 4:
            self.thisptr.refbox = x
        else:
            raise Exception('Reference box should be an iterable of size 4')

    @property
    def bbox(self):
        return self.thisptr.bbox

    @bbox.setter
    def bbox(self, x):
        if len(x) == 4:
            self.thisptr.bbox = x
        else:
            raise Exception('Bounding box should be an iterable of size 4')

    @property
    def demErrorFlag(self):
        return self.thisptr.estimateDEMError

    @demErrorFlag.setter
    def demErrorFlag(self, x):
        self.thisptr.estimateDEMError = x

    @property
    def startingRange(self):
        return self.thisptr.startingRange

    @startingRange.setter
    def startingRange(self, x):
        self.thisptr.startingRange = x

    @property
    def rangeSpacing(self):
        return self.thisptr.rangeSpacing

    @rangeSpacing.setter
    def rangeSpacing(self,x ):
        self.thisptr.rangeSpacing = x

    @property
    def wavelength(self):
        return self.thisptr.wavelength

    @wavelength.setter
    def wavelength(self, x):
        self.thisptr.wavelength=x

    @property
    def incidenceAngleFile(self):
        return self.thisptr.incAngleFile.decode('utf-8')

    @incidenceAngleFile.setter
    def incidenceAngleFile(self, x):
        self.thisptr.incAngleFile = x.encode('utf-8')
