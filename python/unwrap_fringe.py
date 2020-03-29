#!/usr/bin/env python3

# Author: Heresh Fattahi

import os
import argparse
import gdal
import isce
import isceobj

def cmdLineParser():
    '''
    Command line parser.
    '''

    parser = argparse.ArgumentParser(description = 'unwrap estimated wrapped phase',
                formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('-i', '--interferogram_file', type=str, dest='interferogramFile',
            required=True, help='Input interferogram file with complex format')

    parser.add_argument('-c', '--coherence_file', type=str, dest='coherenceFile',
            required=True, help='Input coherence file')

    parser.add_argument('-o', '--unwrap_file', type=str, dest='unwrapFile',
            required=True, help='Output unwrapped file')

    parser.add_argument('-m', '--method', type=str, dest='method',
            default='snaphu', help='unwrapping method: default = snaphu')

    return parser.parse_args()


def unwrap_phass(inps, length, width):

    import Phass

    phass = Phass.pyPhass()
    phass.interferogramFile = inps.interferogramFile
    phass.coherenceFile = inps.coherenceFile
    phass.nr_pixels = width
    phass.nr_lines = length
    phass.outputFile = inps.unwrapFile
    phass.seedFile = inps.unwrapFile + ".seed"
    phass.flowFile = inps.unwrapFile + ".flow"
    phass.visitPatchFile = inps.unwrapFile + ".reg"

    phass.qthresh = 0.1

    phass.unwrap()

    write_xml(phass.outputFile, width, length, 1 , "FLOAT", "BIL")

def unwrap_snaphu(inps, length, width):
    import isce
    import isceobj
    from contrib.Snaphu.Snaphu import Snaphu
   
    altitude = 800000.0
    earthRadius = 6371000.0
    wavelength = 0.056

    snp = Snaphu()
    snp.setInitOnly(False)
    snp.setInput(inps.interferogramFile)
    snp.setOutput(inps.unwrapFile)
    snp.setWidth(width)
    snp.setCostMode('DEFO')
    snp.setEarthRadius(earthRadius)
    snp.setWavelength(wavelength)
    snp.setAltitude(altitude)
    snp.setCorrfile(inps.coherenceFile)
    snp.setInitMethod('MST')
   # snp.setCorrLooks(corrLooks)
    snp.setMaxComponents(100)
    snp.setDefoMaxCycles(2.0)
    snp.setRangeLooks(1)
    snp.setAzimuthLooks(1)
    snp.setCorFileFormat('FLOAT_DATA')
    snp.prepare()
    snp.unwrap()

    write_xml(inps.unwrapFile, width, length, 2 , "FLOAT", "BIL")

def write_xml(fileName,width,length,bands,dataType,scheme):

    img = isceobj.createImage()
    img.setFilename(fileName)
    img.setWidth(width)
    img.setLength(length)
    img.setAccessMode('READ')
    img.bands = bands
    img.dataType = dataType
    img.scheme = scheme
    img.renderHdr()
    img.renderVRT()

    return None

def getSize(data):

    ds = gdal.Open(data, gdal.GA_ReadOnly)
    length = ds.RasterYSize
    width = ds.RasterXSize
    ds = None
    return length, width


if __name__ == '__main__':
    '''
    Main driver.
    '''
    #*************************************************************#
    # read the input options and unwrap
    inps = cmdLineParser()
    length, width = getSize(inps.interferogramFile)
    print("length, width: ", length , " ", width)

    unwrapDir = os.path.dirname(inps.unwrapFile)
    if not os.path.exists(unwrapDir):
        os.makedirs(unwrapDir)

    if inps.method == "snaphu":
        unwrap_snaphu(inps, length, width)

    elif inps.method == "phass":
        unwrap_phass(inps, length, width)

