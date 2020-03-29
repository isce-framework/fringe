#!/usr/bin/python
# This script performs the KS-Test and is derived from Matlab's kstest2. It returns the H value, P value and the KS-test statistic.
# As input it requires two unsorted vectors, the alpha value and the type of test i.e. smaller, larger or unequal.
# The script is invoked with the following $ python pyKstest.py <vec1> <vec2> <alpha>, <test type: smaller, larger, unequal>

import os
import time
import random
import sys
import numpy as np
from collections import OrderedDict

#    # Test input, uncomment the lines below
#x1 = [1, 1, 2, 2, 1, 3, 2]
#x2 = [4, 5, 7, 2, 5, 2, 4]
##
#alpha = 0.05
#tail = "unequal"

def main(x1, x2, alpha=0.05, tail="unequal"):
#def main():

#    if len(sys.argv) != 5:
#        sys.stderr.write('Usage: ./pyKstest.py <vec1> <vec2> <alpha>, <test type: smaller, larger, unequal>\n')
#        sys.exit(1)
#
#    x1 = float(sys.argv[1])
#    x2 = float(sys.argv[2])
#    alpha = float(sys.argv[3])
#    tail = sys.argv[4]
#    print len(x1), len(x2)
#    print x1 #, x2

    binEdges = np.hstack([0, np.sort(np.concatenate([x1,x2])), 100000000])

    binCounts1 = np.histogram(x1, binEdges)[0]
    binCounts2 = np.histogram(x2, binEdges)[0]

    sampleCDF1 = np.cumsum(binCounts1, dtype=float)/np.sum(binCounts1)
    sampleCDF2 = np.cumsum(binCounts2, dtype=float)/np.sum(binCounts2)

    #


    if tail == "unequal":
        deltaCDF = np.abs(np.array(sampleCDF1) - np.array(sampleCDF2))
    elif tail == "smaller":
        deltaCDF = np.array(sampleCDF2) - np.array(sampleCDF1)
    elif tail == "larger":
        deltaCDF = np.array(sampleCDF1) - np.array(sampleCDF2)


    KSstatistic = np.max(deltaCDF)

    n1         = len(x1)
    n2         = len(x2)

    n         = float(n1 * n2) / (n1 + n2)
    lambd    = np.max((np.sqrt(n) + 0.12 + 0.11/np.sqrt(n)) * KSstatistic, 0)
    #print lambd

    if tail != "unequal":
        pValue = np.exp(-2 * lambd * lambd)
    else:
    #    j     =    np.transpose(np.matrix(np.linspace(0,100,num=101)))
        j     =    np.linspace(1,101,num=101)
        pValue = 2 * sum((np.power(-1, j-1) * np.exp(-2 * lambd * lambd * np.power(j, 2)))) # multiple by 2 in MATLABss
        pValue = min(max(pValue,0), 1)

    # assign the H values after calculation is done
    if alpha >= pValue:
        H = 1
    else:
        H = 0

#    print H, pValue, KSstatistic
    return H, pValue, KSstatistic

if __name__ == '__main__':

    files = ['Rayleigh1.bin', 'Rayleigh2.bin',
             'Rayleigh3.bin', 'Exponential.bin']

    for infile in files:
        data = np.fromfile(infile, dtype=np.float32)
        nrows = data.size // 2
        data = data.reshape((2, nrows))

        H, pvalue, kstat = main(data[0,:], data[1,:])

        print('File: {0} P-value: {1}  K-stat: {2}'.format(infile, pvalue, kstat))
        print('Hypothesis: {0}'.format(H))
