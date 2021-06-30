#!/usr/bin/env python3

# Author: Heresh Fattahi

import os
import time
import argparse
from osgeo import gdal
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

    parser.add_argument('-x', '--xml_file', type=str, dest='xmlFile',
            required=False, help='path of reference xml file for unwrapping with snaphu')

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

#Adapted code from unwrap.py and s1a_isce_utils.py in topsStack
def extractInfo(inps):

    '''
    Extract required information from pickle file.
    '''
    from isceobj.Planet.Planet import Planet
    from isceobj.Util.geo.ellipsoid import Ellipsoid
    from iscesys.Component.ProductManager import ProductManager as PM

    pm = PM()
    #pm.configure
    frame = pm.loadProduct(inps.xmlFile)

    burst = frame.bursts[0]
    planet = Planet(pname='Earth')
    elp = Ellipsoid(planet.ellipsoid.a, planet.ellipsoid.e2, 'WGS84')

    data = {}
    data['wavelength'] = burst.radarWavelength

    tstart = frame.bursts[0].sensingStart
    #tend   = frame.bursts[-1].sensingStop
    #tmid = tstart + 0.5*(tend - tstart)
    tmid = tstart


    orbit = burst.orbit
    peg = orbit.interpolateOrbit(tmid, method='hermite')

    refElp = Planet(pname='Earth').ellipsoid
    llh = refElp.xyz_to_llh(peg.getPosition())
    hdg = orbit.getENUHeading(tmid)
    refElp.setSCH(llh[0], llh[1], hdg)

    earthRadius = refElp.pegRadCur

    altitude   = llh[2]

    #sv = burst.orbit.interpolateOrbit(tmid, method='hermite')
    #pos = sv.getPosition()
    #llh = elp.ECEF(pos[0], pos[1], pos[2]).llh()

    data['altitude'] = altitude  #llh.hgt

    #hdg = burst.orbit.getHeading()
    data['earthRadius'] = earthRadius  #elp.local_radius_of_curvature(llh.lat, hdg)
    return data

def unwrap_snaphu(inps, length, width, metadata):
    from contrib.Snaphu.Snaphu import Snaphu

    if metadata is None:
        altitude = 800000.0
        earthRadius = 6371000.0
        wavelength = 0.056
    else:
        altitude = metadata['altitude']
        earthRadius = metadata['earthRadius']
        wavelength = metadata['wavelength']

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
    write_xml(inps.unwrapFile+'.conncomp', width, length, 1, "BYTE", "BIP")


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
    start_time = time.time()

    #*************************************************************#
    # read the input options and unwrap
    inps = cmdLineParser()
    length, width = getSize(inps.interferogramFile)
    print("length, width: ", length , " ", width)

    unwrapDir = os.path.dirname(inps.unwrapFile)
    if not os.path.exists(unwrapDir):
        os.makedirs(unwrapDir)

    if inps.method == "snaphu":
        if inps.xmlFile is not None:
            metadata = extractInfo(inps)
        else:
            metadata = None
        unwrap_snaphu(inps, length, width, metadata)

    elif inps.method == "phass":
        unwrap_phass(inps, length, width)

    m, s = divmod(time.time()-start_time, 60)
    print('time used: {:02.0f} mins {:02.1f} secs.\n'.format(m, s))
