#!/usr/bin/env python3

import os
import argparse

def cmdLineParser():
    '''
    Command line parser.
    '''

    parser = argparse.ArgumentParser(description = 'Create neighborhood mask and count map using KS statistics for stack of coregistered SLCs',
                formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('-i', '--input', type=str, dest='inputDS',
            required=True, help='Input GDAL SLC stack VRT')
    parser.add_argument('-o', '--output', type=str, dest='outputDS',
            required=True, help='Output neighborhood weights mask')
    parser.add_argument('-c', '--count', type=str, dest='countDS',
            required=True, help='Output count dataset')
    parser.add_argument('-m', '--mask', type=str, dest='maskDS',
            default='', help='Optional mask layer to speed up computation')
    parser.add_argument('-l', '--linesperblock', type=int, dest='linesPerBlock',
            default=64, help='Quantum for block of lines')
    parser.add_argument('-r', '--ram', type=int, dest='memorySize',
            default=256, help='Memory in Mb to use')
    parser.add_argument('-x', '--xhalf', type=int, dest='halfWindowX',
            default=5, help='Half window size (range)')
    parser.add_argument('-y', '--yhalf', type=int, dest='halfWindowY',
            default=5, help='Half window size (azimuth)')
    parser.add_argument('-p', '--prob', type=float, dest='pValue',
            default=0.05, help='Minimum p-value for labeling neighbors.')
    parser.add_argument('-s', '--stat', type=str, dest='method',
            default='KS2', help='Statistical test to use - KS2 or AD2')
    parser.add_argument('--nogpu', dest='noGPU', action='store_true',
            default=False, help='Explicitly do not use GPU')

    return parser.parse_args()


def runNmap(inps):
    '''
    Actually run Nmap.
    '''

    import nmaplib

    aa = nmaplib.Nmap()
    
    ###Explicit wiring. Can be automated later.
    aa.inputDS = inps.inputDS
    aa.weightsDS = inps.outputDS
    aa.countDS = inps.countDS
    if (inps.maskDS):
        aa.maskDS = inps.maskDS

    aa.blocksize = inps.linesPerBlock
    aa.memsize = inps.memorySize
    aa.halfWindowX = inps.halfWindowX
    aa.halfWindowY = inps.halfWindowY
    aa.minimumProbability = inps.pValue
    aa.method = inps.method.upper()
    aa.noGPU = inps.noGPU
    aa.run()


if __name__ == '__main__':
    '''
    Main driver.
    '''

    inps = cmdLineParser()

    outDir = os.path.abspath(os.path.dirname(inps.outputDS))
    if not os.path.exists(outDir):
        os.makedirs(outDir)

    runNmap(inps)

