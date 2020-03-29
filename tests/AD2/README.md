# Anderson-Darling 2-sample test

## Scipy implementation

scipy.stats has an anderson\_ksamp method that was used for comparison. This implements the "AakN2" version of the stats.


## R implementation

We used "ad\_test" method in "kSamples" package for comparison. Turns out that this is the same code that was used in CERN's ROOT package almost verbatim. So feel more confident about this.

This script was used to test the python implementation. Since, I do not have access to MATLAB licenses I will assume that this python implementation is equivalent to the MATLAB implementation


## CPP implementation

Implementation adapted from CERN's ROOT Data Analysis Framework [https://root.cern.ch](https://root.cern.ch)

All credit for this module goes to developers  contributors of the ROOT package and R software. 
Source: ROOT/math/mathcore/src/GoFTest.cxx:AndersonDarling2SampleTest

In our implementation, we optimize performance assuming all the samples are unique. Reduces time spent on sorts or bin counting. AD2 test is in general slower than KS2 test as the sorted arrays are again sorted into a single array for computing stats.


## Test data generation

Generated using makeData.py. The true p-values etc are not known. This is meant to be cross comparison between Scipy, R and C++ implementations. 4 different pairs of samples are used. 


## Conclusion


1. There appears to be numerous definitions of the AD2 test - one with AkN2 stat and one with AakN2 stat. We have implemented the first one.

2. R, Scipy and ROOT implementations differ in the implementation of the probability function.

3. Our implementation use ROOT package's P-value function which is not too different but probably a slightly different approximation. 

4. The K-statistic computed by all the different implementations are the same.

5. Since, in this work we will work with finite number of SAR images in a stack we will stick to the ROOT implementation. This does not seem to be too different from R or Scipy for most cases for hypothesis testing.
