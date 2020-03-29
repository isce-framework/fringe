#!/usr/bin/env python3

import os
import argparse

def cmdLineParser():
    '''
    Command line parser.
    '''

    parser = argparse.ArgumentParser(description = 'Compute amplitude calibration using magnitude of SLCs',
                formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('-i', '--input', type=str, dest='inputDS',
            required=True, help='Input GDAL SLC stack VRT')
    parser.add_argument('-o', '--output', type=str, dest='outputDS',
            required=True, help='Output GDAL SLC stack VRT')
    parser.add_argument('-m', '--mask', type=str, dest='maskDS',
            default='', help='Mask DS to use for masking out water / other areas for estimating calibration constant.')
    parser.add_argument('-d', '--default', type=float, dest='defaultValue',
            default=1.0, help='Default value to assign to mask if a constant cannot be computed')
    parser.add_argument('-l', '--linesperblock', type=int, dest='linesPerBlock',
            default=64, help='Quantum for block of lines')
    parser.add_argument('-r', '--ram', type=int, dest='memorySize',
            default=256, help='Memory in Mb to use')
    parser.add_argument('-s', '--sqrt', action='store_true', default=False,
            help='Apply sqrt explicitly to the amplitudes')

    return parser.parse_args()


def runCalamp(indict):
    '''
    Actually run calamp.
    '''

    import calamplib

    aa = calamplib.Calamp()
    
    ###Explicit wiring. Can be automated later.
    aa.inputDS = inps.inputDS
    aa.outputDS = inps.outputDS

    if (inps.maskDS):
        aa.maskDS = inps.maskDS

    aa.defaultValue = inps.defaultValue
    aa.blocksize = inps.linesPerBlock
    aa.memsize = inps.memorySize
    aa.applySqrt = inps.sqrt

    aa.run()


if __name__ == '__main__':
    '''
    Main driver.
    '''

    inps = cmdLineParser()

    runCalamp(inps)

