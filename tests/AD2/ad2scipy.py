#!/usr/bin/python

import os
import time
import random
import sys
import numpy as np
from scipy.stats import anderson_ksamp

if __name__ == '__main__':

    files = ['Rayleigh1.bin', 'Rayleigh2.bin',
             'Rayleigh3.bin', 'Exponential.bin']

    for infile in files:
        data = np.fromfile(infile, dtype=np.float32)
        nrows = data.size // 2
        data = data.reshape((2, nrows))

        stat, crit, pvalue = anderson_ksamp([data[0],data[1]])

        print('File: {0} P-value: {1}  AD-stat: {2}'.format(infile, pvalue, stat))
