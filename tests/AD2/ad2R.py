#!/usr/bin/python

import os
import time
import random
import sys
import numpy as np
from collections import OrderedDict
import rpy2.robjects as R
from rpy2.robjects.packages import importr

if __name__ == '__main__':



    stats = importr('kSamples')

    files = ['Rayleigh1.bin', 'Rayleigh2.bin',
             'Rayleigh3.bin', 'Exponential.bin']

    for infile in files:
        data = np.fromfile(infile, dtype=np.float32)
        nrows = data.size // 2
        data = data.reshape((2, nrows))

        x1 = R.FloatVector([float(x) for x in data[0]])
        x2 = R.FloatVector([float(x) for x in data[1]])

        v = stats.ad_test(x1,x2,method='exact')
        mat = np.array(v[6])
        kstat = mat[0][1]
        pvalue = mat[0][3]

        print('File: {0} P-value: {1}  K-stat: {2}'.format(infile, pvalue, kstat))
