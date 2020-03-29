#!/usr/bin/env python3

import numpy as np 
import matplotlib.pyplot as plt
def wrap_phase(ph):
    return ph - 2 * np.pi * np.round(ph / (2*np.pi))

if __name__ == '__main__':
    '''
    Create datasets for testing.
    '''

    ###Example 1
    wvl = 0.06
    rng = 850000.
    inc = 34.0
    delz = 6.0 
    nifg = 30
    bperp = np.random.randn(nifg) * 1000.
    ph =  4 * np.pi * bperp * delz /(wvl * rng * np.sin( np.radians(inc))) +  np.random.random(nifg)
    ph = wrap_phase(ph)

    fid = open('test1.dat', 'wb')
    np.float32(wvl).tofile(fid)
    np.float32(rng).tofile(fid)
    np.float32(inc).tofile(fid)
    np.float32(delz).tofile(fid)

    bperp.astype(np.float32).tofile(fid)
    ph.astype(np.float32).tofile(fid)


    ###Example 2
    wvl = 0.06
    rng = 920000.
    inc = 40.0
    delz = 20.0 
    nifg = 35
    bperp = np.random.randn(nifg) * 500.
    ph =  4 * np.pi * bperp * delz /(wvl * rng * np.sin( np.radians(inc))) +  np.random.random(nifg)
    ph = wrap_phase(ph)

    fid = open('test2.dat', 'wb')
    np.float32(wvl).tofile(fid)
    np.float32(rng).tofile(fid)
    np.float32(inc).tofile(fid)
    np.float32(delz).tofile(fid)

    bperp.astype(np.float32).tofile(fid)
    ph.astype(np.float32).tofile(fid)

    ###Example 3
    wvl = 0.03
    rng = 788000.
    inc = 23.0
    delz = -10.0 
    nifg = 25
    bperp = np.random.randn(nifg) * 700.
    ph =  4 * np.pi * bperp * delz /(wvl * rng * np.sin( np.radians(inc))) +  np.random.random(nifg)
    ph = wrap_phase(ph)

    fid = open('test3.dat', 'wb')
    np.float32(wvl).tofile(fid)
    np.float32(rng).tofile(fid)
    np.float32(inc).tofile(fid)
    np.float32(delz).tofile(fid)

    bperp.astype(np.float32).tofile(fid)
    ph.astype(np.float32).tofile(fid)

