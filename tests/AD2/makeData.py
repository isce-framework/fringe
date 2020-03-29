#!/usr/bin/env python3

import numpy as np 

if __name__ == '__main__':
    '''
    Create datasets for testing.
    '''

    ###Example 1: Rayleigh with same scale
    x1 = np.random.rayleigh(scale=1.0, size=30).astype(np.float32)
    x2 = np.random.rayleigh(scale=1.0, size=30).astype(np.float32)

    x1.sort()
    x2.sort()

    with open('Rayleigh1.bin', 'wb') as fid:
        x1.tofile(fid)
        x2.tofile(fid)


    ###Example 2: Rayleigh with slightly different scale
    x1 = 20*np.random.rayleigh(scale=1.0, size=50).astype(np.float32)
    x2 = 20*np.random.rayleigh(scale=1.5, size=50).astype(np.float32)

    x1.sort()
    x2.sort()

    with open('Rayleigh2.bin', 'wb') as fid:
        x1.tofile(fid)
        x2.tofile(fid)

   
    ###Example 3: Rayleigh with different scales
    x1 = 20*np.random.rayleigh(scale=1.0, size=40).astype(np.float32)
    x2 = 30*np.random.rayleigh(scale=2.5, size=40).astype(np.float32)

    x1.sort()
    x2.sort()

    with open('Rayleigh3.bin', 'wb') as fid:
        x1.tofile(fid)
        x2.tofile(fid)


    ####Example 4: Same exponential distribution
    x1 = 10*np.random.exponential(scale=1.5, size=25).astype(np.float32)
    x2 = 10*np.random.exponential(scale=1.5, size=25).astype(np.float32)

    x1.sort()
    x2.sort()

    with open('Exponential.bin', 'wb') as fid:
        x1.tofile(fid)
        x2.tofile(fid)

