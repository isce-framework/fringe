#!/usr/bin/env python3

import argparse
import numpy as np
from osgeo import gdal
import pandas 
import os

def cmdLineParse():
    '''
    Command line parser.
    '''

    parser = argparse.ArgumentParser(description='SBAS analysis')
    parser.add_argument('--input', dest='infile', type=str,
                        required=True, help='Input csv file with interferograms that Pandas can read')
    parser.add_argument('--baseline', dest='bperpfile', type=str,
                        default=None, help='Input csv file with interferograms that Pandas can read')
    parser.add_argument('--outdir', dest='outdir', type=str,
                        required=True, help='Output directory')
    parser.add_argument('--ref', dest='refbox', type=int,
                        nargs=4, required=True,
                        help='Reference box in xoff, yoff, nx, ny format')
    parser.add_argument('--box', dest='bbox', type=int,
                        nargs=4, default=[0, 0, 0, 0],
                        help='Portion of image to process in xoff, yoff, nx, ny format')
    parser.add_argument('--vrt', dest='vrtdir', type=str,
                        default='VRT', help='VRT pointing to original data')

    parser.add_argument('--thresh', dest='thresh', type=float,
                        default=0.2, help='Quality threshold')

    parser.add_argument('--scale', dest='scale', type=float,
                        default=1.0, help='Scale factor')
    parser.add_argument('--date', dest='refdate', type=str,
                        default=None, help='Reference date')
    parser.add_argument('--wvl', dest='wavelength', type=float,
                        default=0.0, help='Radar wavelength')
    parser.add_argument('--startingrange', dest='startingrange', type=float,
                        default=0.0, help='Starting range')
    parser.add_argument('--rangespacing', dest='rangespacing', type=float,
                        default=0.0, help='Range spacing')
    parser.add_argument('--incangle', dest='incangle', type=str,
                        default='', help='Incidence angle file')

    parser.add_argument('--demerr', dest='estimateDEMError', action='store_true',
                        default=False, help='Estimate DEM Error flag')
    inps = parser.parse_args()

    return inps

def loadInput(csvfile, ifg=False):
    '''
    Load csv file using Pandas
    '''
    import pandas

    df = pandas.read_csv(csvfile)
    
    ##Change all columns to lower case
    df.columns = [x.lower() for x in df.columns]
    if ifg:
        df.set_index(['master date', 'slave date'])

    return df

def generateVRT(dest, src):
    '''
    Create VRT using Python API.
    '''

    ds = gdal.Open(src, gdal.GA_ReadOnly)
    opts = gdal.TranslateOptions(format='VRT', outputType=gdal.GDT_Float32,
                                 bandList=[ds.RasterCount])
                                 
    gdal.Translate(dest, ds, options=opts)
   
    ds = None

    return dest


def createSBAS(df, inps, sardf=None):
    '''
    Create SBAS object from Pandas dataframe
    '''
    
    def setIfValid(src, default=None):
        if not pandas.isnull(src):
            return src
        elif default is not None:
            return default
    
    def cannotBeInvalid(src):
        if pandas.isnull(src):
            raise Exception('Invalid NULL entry for mandatory data field')
        return src

    import sbaslib

    sbasObj = sbaslib.SBASWrapper()

    nPairs = len(df.index)
    sbasObj.setNumberOfPairs( nPairs)
    print('Number of pairs: ', nPairs)

    sbasObj.referenceBox = inps.refbox
    sbasObj.bbox = inps.bbox

    ##Setup output directory
    if not os.path.exists(inps.outdir):
        os.makedirs(inps.outdir)

    else:
        raise Exception('Output Dir {0} already exists. Exitting ...'.format(inps.outdir))

    sbasObj.outputDir = inps.outdir

    ###Check if data frame has columns
    hasThreshold = 'threshold' in df.columns
    hasScale = 'scale' in df.columns

    if not os.path.exists(inps.vrtdir):
        os.makedirs(inps.vrtdir)

    ###Start wiring 
    for ii in range(nPairs):
        pairIndex = df.index[ii]
        pairObj = sbasObj.getPair(pairIndex)
        rec = df.loc[pairIndex]

        pairObj.masterDate = str(cannotBeInvalid(rec['master date']))
        pairObj.slaveDate = str(cannotBeInvalid(rec['slave date']))
        pairObj.ifgName = generateVRT( os.path.join(inps.vrtdir, '{0}_{1}.unw.vrt'.format(pairObj.masterDate, pairObj.slaveDate)), cannotBeInvalid( rec['unwrapped interferogram']))

        pairObj.cohName = generateVRT(os.path.join(inps.vrtdir, '{0}_{1}.coh.vrt'.format(pairObj.masterDate, pairObj.slaveDate)), cannotBeInvalid(rec['coherence'])) 
        if hasThreshold:
           pairObj.threshold = setIfValid(rec['threshold'],
                        default=inps.thresh)

        if hasScale:
            pairObj.scale = setIfValid(rec['scale'],
                        default=inps.scale)


    ###Check inputs and set up as needed
    if inps.refdate is not None:
        sbasObj.referenceDate = inps.refdate
    

    if inps.estimateDEMError:
        sbasObj.demErrorFlag = True

        if any([x==0. for x in [inps.startingrange, inps.rangespacing, inps.wavelength]]):
            raise Exception('DEM Error Estimation requested but range inputs not provided')

        if inps.incangle in [None, '']:
            raise Exception('DEM Error Estimation requested but incidence angle file not provided')

        sbasObj.startingRange = inps.startingrange
        sbasObj.rangeSpacing = inps.rangespacing
        sbasObj.wavelength = inps.wavelength
        sbasObj.incidenceAngleFile = inps.incangle

        if sardf is None:
            raise Exception('No SAR baselines provided for dem error estimation')

        sbasObj.setNumberOfScenes(len(sardf.index))
        for ii in range(len(sardf.index)):
            sceneIndex = sardf.index[ii]
            sceneObj = sbasObj.getScene(ii)
            rec = sardf.loc[sceneIndex]

            sceneObj.date = str(rec['date'])
            sceneObj.baselineFile = rec['baseline']

    return sbasObj

if __name__ == '__main__':
    '''
    Main driver.
    '''
    #Parse command line
    inps = cmdLineParse()

    #Read in interferogram data
    df = loadInput(inps.infile, ifg=True)

    if inps.bperpfile not in ['', None]:
        sardf = loadInput(inps.bperpfile)
    else:
        sardf = None

    #Create SBAS object
    sbasObj = createSBAS(df, inps, sardf=sardf)

    #Process
    sbasObj.process()
