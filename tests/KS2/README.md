# Kolmogorov-Smirnov 2-sample test

## Matlab implementation

Matlab script derived from: 
[https://github.com/ICEACE/MATLAB/blob/master/kstest2.m](https://github.com/ICEACE/MATLAB/blob/master/kstest2.m)

This script was used to test the matlab implementation.

## Python implementation

An equivalent Python script was found here:
[https://github.com/Zahaib/Python-Kstest2/blob/master/pyKstest.py](https://github.com/Zahaib/Python-Kstest2/blob/master/pyKstest.py)

This script was used to test the python implementation. Since, I do not have access to MATLAB licenses I will assume that this python implementation is equivalent to the MATLAB implementation


## CPP implementation

Implementation adapted from CERN's ROOT Data Analysis Framework [https://root.cern.ch](https://root.cern.ch)

All credit for this module goes to developers  contributors of the ROOT package. 
Source: ROOT/math/mathcore/src/TMath.cxx:KolmogorovTest


## Test data generation

Generated using makeData.py. The true p-values etc are not known. This is meant to be cross comparison between MATLAB, Python and C++ implementations. 4 different pairs of samples are used. 


## Conclusion

1. MATLAB appears to have implemented an asymptotic P-value function that is good for large sample sizes using "lambda" variable.

2. R appears to have implemented a different P-value function based on the actual K-statistic.

3. Our implementation use ROOT package's P-value function which is not too different but probably a slightly different approximation. 

4. The K-statistic computed by all the different implementations are the same.

5. Since, in this work we will work with finite number of SAR images in a stack we will stick to the ROOT implementation. This does not seem to be too different from R or MATLAB for most cases for hypothesis testing.


The following thread on stack overflow also appears to have reached the same conclusion regarding P-values and asymptotic P-values:
[https://stackoverflow.com/questions/38933425/poorly-implemented-two-sample-kolmogorov-smirnov-test-kstest2-in-matlab](https://stackoverflow.com/questions/38933425/poorly-implemented-two-sample-kolmogorov-smirnov-test-kstest2-in-matlab)

