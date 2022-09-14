#!/usr/bin/env python3

import os
import argparse

def cmdLineParser():
    '''
    Command line parser.
    '''

    parser = argparse.ArgumentParser(description = 'Perform MLE-based phase-linking on a stack of coregistered SLCs',
                formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('-i', '--input', type=str, dest='inputDS',
            required=True, help='Input GDAL SLC stack VRT')
    parser.add_argument('-w', '--wts', type=str, dest='wtsDS',
            required=True, help='Input weights dataset')
    parser.add_argument('-o', '--output', type=str, dest='outputFolder',
            required=True, help='Output folder with phase-linked SLCs')
    parser.add_argument('-l', '--linesperblock', type=int, dest='linesPerBlock',
            default=64, help='Quantum for block of lines')
    parser.add_argument('-r', '--ram', type=int, dest='memorySize',
            default=2048, help='Memory in Mb to use')
    parser.add_argument('-x', '--xhalf', type=int, dest='halfWindowX',
            default=5, help='Half window size (range)')
    parser.add_argument('-y', '--yhalf', type=int, dest='halfWindowY',
            default=5, help='Half window size (azimuth)')
    parser.add_argument('-n', '--minneigh', type=int, dest='minNeighbors',
            default=5, help='Minimum number of neighbors for computation')
    parser.add_argument('-m', '--method', type=str, dest='method',
            default='MLE', help='Decomposition method to use - MLE / EVD / STBAS')
    parser.add_argument('-b', '--bandwidth', type=int, dest='bandWidth',
            default=-1, help='Diagonal bandwidth for STBAS')

    return parser.parse_args()


def runEvd(inps):
    '''
    Actually run Evd.
    '''

    import os
    os.environ['OPENBLAS_NUM_THREADS'] = "1"

    import phase_linklib
    aa = phase_linklib.Phaselink()
    
    ###Explicit wiring. Can be automated later.
    aa.inputDS = inps.inputDS
    aa.weightsDS = inps.wtsDS
    aa.outputFolder = inps.outputFolder

    aa.outputCompressedSlcFolder = aa.outputFolder
    aa.compSlc = "compslc.bin"

    aa.blocksize = inps.linesPerBlock
    aa.memsize = inps.memorySize
    aa.halfWindowX = inps.halfWindowX
    aa.halfWindowY = inps.halfWindowY
    aa.minimumNeighbors = inps.minNeighbors
    
    ##Set up method and bandwidth
    aa.method = inps.method
    aa.bandWidth = inps.bandWidth

    aa.run()


if __name__ == '__main__':
    '''
    Main driver.
    '''

    inps = cmdLineParser()

    runEvd(inps)

