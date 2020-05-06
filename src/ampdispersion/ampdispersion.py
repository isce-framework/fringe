#!/usr/bin/env python3

import os
import argparse

def cmdLineParser():
    '''
    Command line parser.
    '''

    parser = argparse.ArgumentParser(description = 'Compute amplitude dispersion and mean amplitude for stack of coregistered SLCs',
                formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('-i', '--input', type=str, dest='inputDS',
            required=True, help='Input GDAL SLC stack VRT')
    parser.add_argument('-o', '--output', type=str, dest='outputDS',
            required=True, help='Output amplitude dispersion dataset')
    parser.add_argument('-m', '--mean', type=str, dest='meanampDS',
            default='', help='Output mean amplitude')
    parser.add_argument('-l', '--linesperblock', type=int, dest='linesPerBlock',
            default=64, help='Quantum for block of lines')
    parser.add_argument('-r', '--ram', type=int, dest='memorySize',
            default=256, help='Memory in Mb to use')
    parser.add_argument('-b', '--band', type=int, dest='refBand',
            default=1, help='Reference band to use for relative normalization')

    return parser.parse_args()


def runAmpdispersion(indict):
    '''
    Actually run ampdispersion.
    '''

    import ampdispersionlib

    aa = ampdispersionlib.Ampdispersion()
    
    ###Explicit wiring. Can be automated later.
    aa.inputDS = inps.inputDS
    aa.outputDS = inps.outputDS
    aa.meanampDS = inps.meanampDS

    aa.blocksize = inps.linesPerBlock
    aa.memsize = inps.memorySize
    aa.refband = inps.refBand

    aa.run()

    return


if __name__ == '__main__':
    '''
    Main driver.
    '''

    inps = cmdLineParser()

    outDir = os.path.abspath(os.path.dirname(inps.outputDS))
    if not os.path.exists(outDir):
        os.makedirs(outDir)

    runAmpdispersion(inps)

    # create xml file if missing
    for fname in [inps.outputDS, inps.meanampDS]:
        if not os.path.isfile(fname+'.xml'):
            cmd = 'gdal2isce_xml.py -i {}'.format(fname)
            print(cmd)
            os.system(cmd)


