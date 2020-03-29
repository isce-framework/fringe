#!/usr/bin/python

import os
import time
import random
import sys
import numpy as np
from scipy.stats import ks_2samp

if __name__ == '__main__':

    files = ['Rayleigh1.bin', 'Rayleigh2.bin',
             'Rayleigh3.bin', 'Exponential.bin']

    for infile in files:
        data = np.fromfile(infile, dtype=np.float32)
        nrows = data.size // 2
        data = data.reshape((2, nrows))

        kstat, pvalue = ks_2samp(data[0],data[1])

        print('File: {0} P-value: {1}  K-stat: {2}'.format(infile, pvalue, kstat))
