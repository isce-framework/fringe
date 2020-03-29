#!/usr/bin/env python3

import os
import argparse

def cmdLineParser():
    '''
    Command line parser.
    '''

    parser = argparse.ArgumentParser(description = 'Despeckle SLC amplitude / single look interferogram',
                formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('-i', '--input', type=str, dest='inputDS',
            required=True, help='Input GDAL SLC stack VRT')
    parser.add_argument('-o', '--output', type=str, dest='outputDS',
            required=True, help='Output despeckled dataset')
    parser.add_argument('-w', '--wts', type=str, dest='wtsDS',
            default='', help='Optional mask layer to speed up computation')
    parser.add_argument('-l', '--linesperblock', type=int, dest='linesPerBlock',
            default=64, help='Quantum for block of lines')
    parser.add_argument('-r', '--ram', type=int, dest='memorySize',
            default=512, help='Memory in Mb to use')
    parser.add_argument('-x', '--xhalf', type=int, dest='halfWindowX',
            default=5, help='Half window size (range)')
    parser.add_argument('-y', '--yhalf', type=int, dest='halfWindowY',
            default=5, help='Half window size (azimuth)')
    parser.add_argument('-b', '--band', type=int, dest='bands', nargs='*', default=[],
            help='Single band for SLC despeckling and 2 bands for interferogram')
    parser.add_argument('-c', '--corr', action='store_true', dest='cohFlag', default=False,
            help='Replace amplitude for despecked interferogram with coherence')

    return parser.parse_args()


def runDespeck(inps):
    '''
    Actually run despeck.
    '''

    import despecklib

    aa = despecklib.Despeck()
    
    ###Explicit wiring. Can be automated later.
    aa.inputDS = inps.inputDS
    aa.weightsDS = inps.wtsDS
    aa.outputDS = inps.outputDS
    aa.blocksize = inps.linesPerBlock
    aa.memsize = inps.memorySize
    aa.halfWindowX = inps.halfWindowX
    aa.halfWindowY = inps.halfWindowY

    if len(inps.bands) == 1:
        aa.band1 = inps.bands[0]
        if inps.cohFlag:
            raise Exception('User requested coherence when requesting despeckling of SLC magnitude')
    elif len(inps.bands) == 2:
        aa.band1 = inps.bands[0]
        aa.band2 = inps.bands[1]
        aa.coherenceFlag = inps.cohFlag
    elif len(inps.bands) != 0:
        raise Exception('Despeck can handle one or two bands. More than two bands provided')

    aa.run()


if __name__ == '__main__':
    '''
    Main driver.
    '''

    inps = cmdLineParser()

    runDespeck(inps)

