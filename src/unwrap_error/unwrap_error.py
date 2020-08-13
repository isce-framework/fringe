#!/usr/bin/env python3
# ----------------------------------------------------------------------------
# Copyright (c) 2018-, California Institute of Technology ("Caltech"). 
# U.S. Government sponsorship acknowledged.
# All rights reserved.
# 
# Author(s): Heresh Fattahi
# ----------------------------------------------------------------------------

import os
import argparse
import datetime
import time
from osgeo import gdal
import numpy as np
import DesignMatrix
import unwraplib as Unwrap

def cmdLineParser():
    '''
    Command line parser.
    '''

    parser = argparse.ArgumentParser(description = 'Fitting a time-series with a mathematical model',
                formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('-r', '--region_dataset', type=str, dest='regionDS',
            required=True, help='Input region dataset (A single VRT, whose bands are pairs)')
    parser.add_argument('-u', '--unwrapped_dataset', type=str, dest='unwDS',
            required=True, help='Input dataset of unwrapped interferograms')
    parser.add_argument('-b', '--blockSize', type=int, dest='blockSize',
            default=128, help='number of lines for each processing block of data')
    parser.add_argument('-o', '--output', type=str, dest='outputDir',
            required=True, help='output directory to store estimated ambiguities')

    return parser.parse_args()   
 
def extractCommonRegions(inps):

    tempMaskFile = "/u/k-data/fattahi/test_unw_err/tempMask.bin"
    tempRegionFile = "tempRegionMap.bin"
    unw = Unwrap.PyUnwrap()
    unw.blocksize = inps.blockSize
    unw.Set_connComp_dataset(inps.regionDS)
    unw.maskDS = tempMaskFile
    unw.ComputeCommonMask()
  
    relabel_regions(tempMaskFile, tempRegionFile)

    return tempRegionFile

def relabel_regions(maskFile, regionFile):
   
    from skimage import measure
    ds = gdal.Open(maskFile, gdal.GA_ReadOnly)
    mask = ds.GetRasterBand(1).ReadAsArray()
    ds = None
    regions  = measure.label(mask, background=0)
   
    region_indices = np.unique(regions)
    print("Total number of regions : ", len(region_indices))

    driver = gdal.GetDriverByName('ENVI')
    ds = driver.Create(regionFile, regions.shape[1], regions.shape[0], 1, 4)
    b = ds.GetRasterBand(1).WriteArray(regions)
    ds = None    
 

def extract_time(inputDS):
    ds = gdal.Open(inputDS, gdal.GA_ReadOnly)
    nBands = ds.RasterCount
    time_array=[]
    for bnd in range(nBands):
        dd = ds.GetRasterBand(bnd+1).GetMetadata("slc")['Date']
        dd = datetime.datetime(*time.strptime(dd, "%Y%m%d")[0:6])
        time_array.append(dd) 

    return np.array(time_array)

def getPairs(stack):
    pairs = []
    pairs_band_dict = {}
    ds = gdal.Open(stack, gdal.GA_ReadOnly)
    nbands = ds.RasterCount
    for bnd in range(nbands):
        dd = ds.GetRasterBand(bnd+1).GetMetadata("")['Dates'] 
        pairs.append(dd)
        pairs_band_dict[dd] = bnd+1
    ds = None
    return pairs, pairs_band_dict

def getG(stack):

    pairs, pairs_band_dict = getPairs(stack)
    dm = DesignMatrix.DesignMatrix()
    dm.pairs = pairs
    dm.configure()
    dm.closure()
    G = dm.getG()
    return pairs, pairs_band_dict, G

def extractData(G, stack, xx, yy, numPixels=10):
    if len(xx)<numPixels:
        numPixels = len(xx)
    ds = gdal.Open(stack, gdal.GA_ReadOnly)
    data = ds.ReadAsArray(int(xx[0]),int(yy[0]),1,1)[:,0,0]
    reference_phase = np.zeros_like(data)
    for b in range(ds.RasterCount):
        reference_phase[b] = ds.GetRasterBand(b+1).GetOffset()

    phase_closure = np.dot(G, data - reference_phase)

    for i in range(1,numPixels):
        data = ds.ReadAsArray(int(xx[i]),int(yy[i]),1,1)[:,0,0]        
        phase_closure = np.hstack((phase_closure, np.dot(G, data - reference_phase)))

    ds = None
    print(phase_closure.shape)
    return phase_closure, numPixels    

def estimateAmbiguities(stack, regionFile):
    # get the design matrix based on closure
    pairs, pairs_band_dict, G = getG(stack)
    C = -2.0*np.pi*G
    print("G", G.shape)
    print("G rank: ", np.linalg.matrix_rank(G))

    # read the region map
    ds = gdal.Open(regionFile, gdal.GA_ReadOnly)
    regions = ds.GetRasterBand(1).ReadAsArray()
    ds = None
    region_indices = np.unique(regions)
    number_regions = len(region_indices) 
    number_pairs = len(pairs)
    print("Total number of regions : ", number_regions)
    print("Total number of pairs : ", number_pairs) 
    # An array to store the estimated interger numbers of jumps for the entire stack
    # each row of this array is the estimated correction vector for a region in the stack
    Ustack = np.zeros((number_regions, number_pairs), dtype = np.int8)

    x,y = np.meshgrid(range(regions.shape[1]),range(regions.shape[0]))    
    for reg in region_indices:
        #reg = 140
        if reg==0:
           continue
        print("region {0} out of {1} ".format(reg, number_regions)) 
        regIdx = regions == reg
        xx = x[regIdx]
        yy = y[regIdx]
        L, numPixels = extractData(G, stack, xx, yy, numPixels=20)
        CC = C.copy()
        for ii in range(1,numPixels):
            CC = np.vstack((CC,C))
        X = np.dot(np.linalg.pinv(CC), L)
        U = np.round(X)
        Ustack[reg,:] = U

    return Ustack, pairs, pairs_band_dict

def write_ambiguity(data, outName):
    driver = gdal.GetDriverByName('GTIFF')
    dst_options = ['COMPRESS=LZW']
    ds = driver.Create(outName , data.shape[1], data.shape[0], 1, gdal.GDT_Byte, dst_options)
    bnd = ds.GetRasterBand(1)
    bnd.WriteArray(data)
    bnd.FlushCache()
    ds = None

def scale_2PI_vrt(tifFile, outName):
    ds = gdal.Open(tifFile, gdal.GA_ReadOnly)
    xsize = ds.RasterXSize
    ysize = ds.RasterYSize
    ds = None

    tmpl = '''<VRTDataset rasterXSize="{xsize}" rasterYSize="{ysize}">
  <VRTRasterBand dataType="Float32" band="1">
  <ComplexSource>
      <SourceFilename relativeToVRT="1">{srcFile}</SourceFilename>
      <SourceBand>1</SourceBand>
      <SourceProperties RasterXSize="{xsize}" RasterYSize="{ysize}" DataType="Byte"/>
      <ScaleOffset>0</ScaleOffset>
      <ScaleRatio>6.28319</ScaleRatio>
    </ComplexSource>
  </VRTRasterBand>
</VRTDataset>'''
    
    with open( outName , 'w') as fid:
        fid.write( tmpl.format(xsize = xsize,
                               ysize = ysize,
                               srcFile = tifFile))
    
def correct_igram_vrt(srcFile, correctionFile, correctedFile):

    ds = gdal.Open(srcFile, gdal.GA_ReadOnly)
    xsize = ds.RasterXSize
    ysize = ds.RasterYSize
    ds = None

    tmpl='''<VRTDataset rasterXSize="{xsize}" rasterYSize="{ysize}">
  <VRTRasterBand dataType="Float32" band="1" subClass="VRTDerivedRasterBand">
    <Description>Sum</Description>
    <PixelFunctionType>sum</PixelFunctionType>
    <SimpleSource>
      <SourceFilename>{sourceFile}</SourceFilename>
      <SourceBand>2</SourceBand>
    </SimpleSource>
    <SimpleSource>
      <SourceFilename>{correctionFile}</SourceFilename>
    </SimpleSource>
  </VRTRasterBand>
</VRTDataset>'''
    with open( correctedFile , 'w') as fid:
        fid.write( tmpl.format(xsize = xsize,
                               ysize = ysize,
                               sourceFile = srcFile,
                               correctionFile = correctionFile))


def assign_pair_region_ambiguity(unwStack, regionStack, commonStackRegionFile, Ustack, pairs, pairs_band_dict, ambiguityDir, correctedDir):

    ds = gdal.Open(commonStackRegionFile, gdal.GA_ReadOnly)
    commonMask = ds.GetRasterBand(1).ReadAsArray()
    ds = None

    ####################################
    # in a given interferogram, each region may overlap with multiple regions in the 
    # common region map. We need to evaluate each interferogram and assign one estimated ambiguity index to each region. 
    regionDS = gdal.Open(regionStack, gdal.GA_ReadOnly)
    unwDS = gdal.Open(unwStack, gdal.GA_ReadOnly)

    for ii, pair in enumerate(pairs):
        pair_region = regionDS.GetRasterBand(pairs_band_dict[pair]).ReadAsArray()
        pair_region_U = np.zeros_like(pair_region)
        pair_region_indices = np.unique(pair_region)
        pair_region_indices = pair_region_indices[pair_region_indices > 0]
        for reg in pair_region_indices:
            # find the common regions between this region and the common stack region file
            # first make a mask for the current region            
            mask = (pair_region==reg).astype(int)

            # multiply the mask by the commonMask file
            mask = commonMask*mask
            
            # multiple regions in the common Mask may overlap with the pairs region map
            # the indices of all regions in the commonMask which overlap with the current pair region
            commonMaskIdx = np.unique(mask)
            commonMaskIdx = commonMaskIdx[commonMaskIdx>0].astype(int)
            U = Ustack[commonMaskIdx, ii]
            pair_region_U[pair_region==reg] = np.round(np.median(U))
        
        # write the estimated U for this pair to file
        outName = os.path.join(ambiguityDir, pair + ".tif")
        correctionFile = os.path.join(ambiguityDir, pair + ".vrt")
        correctedFile = os.path.join(correctedDir, pair + ".vrt")
        write_ambiguity(pair_region_U, outName)  
        # create a vrt file which scales the ambiguity integer number by 2PI
        scale_2PI_vrt(outName, correctionFile)
        
        #get the original unwrapped file
        srcFile = unwDS.GetFileList()[pairs_band_dict[pair]]    
        # use pixel function VRT to add the estimated ambiguity to  the unwrapped file          
        correct_igram_vrt(srcFile, correctionFile, correctedFile) 

    regionDS = None
    unwDS = None

if __name__ == '__main__':
    '''
    Main driver.
    '''
    inps = cmdLineParser()
    inps.regionDS = os.path.abspath(inps.regionDS)
    inps.outputDir = os.path.abspath(inps.outputDir)
    ambiguityDir = os.path.join(inps.outputDir, "ambiguities")
    correctedDir = os.path.join(inps.outputDir, "corrected_pairs")

    if not os.path.exists(ambiguityDir):
        os.makedirs(ambiguityDir)

    if not os.path.exists(correctedDir):
        os.makedirs(correctedDir)    

    regionFile = extractCommonRegions(inps)
    Ustack, pairs, pairs_band_dict = estimateAmbiguities(inps.unwDS,  regionFile)
    assign_pair_region_ambiguity(inps.unwDS, inps.regionDS, regionFile, Ustack, pairs, pairs_band_dict, ambiguityDir, correctedDir)
    #assign_pair_region_ambiguity(inps.regionDS, regionFile, Ustack, pairs, pairs_band_dict, inps.outputDir)
