#!/usr/bin/env python3

# Author: Heresh Fattahi

import os
import glob
import argparse
import numpy as np
from osgeo import gdal
import isce
import isceobj
import time
import datetime

#from Stack import Stack, MiniStack

def cmdLineParser():
    '''
    Command line parser.
    '''

    parser = argparse.ArgumentParser(description = 'Creates VRT files to adjusts mini-stack unwrapped phase series with datum adjustment phase.',
                formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('-s', '--slcDir', type=str, dest='slcDir',
            required=True, help='Input folder which contains vrt files of all slcs')
    
    parser.add_argument('-m', '--miniStackDir', type=str, dest='miniStackDir',
            required=True, help='Input mini stack directory')
 
    parser.add_argument('-d', '--datumDir', type=str, dest='datumDir',
            required=True, help='Input mini stack directory')

    parser.add_argument('-M', '--miniStackSize', type=int, dest='miniStackSize',
            required=True, help='size of each miniStack')

    parser.add_argument('-o', '--outDir', type=str, dest='outDir',
            required=True, help='output directory')

    parser.add_argument('-t', '--timeSeries_vrt_file', type=str, dest='timeSeriesName',
            required=False, help='output vrt file for time-series')

    parser.add_argument('-x', '--reference_x', type=int, dest='ref_x',
            required=False, help='x coordinate of refernce pixel')

    parser.add_argument('-y', '--reference_y', type=int, dest='ref_y',
            required=False, help='y coordinate of refernce pixel')

    parser.add_argument('--unwrapped', action='store_true',
            default=False, help='adjusting unwrapped dataset. Default: False which means adjustment is done on wrapped data')

    return parser.parse_args()

def getDates(slcDir):
    dateList = []

    slcs = glob.glob(os.path.join(slcDir, "*.vrt"))
    for slc in slcs:
        dd = os.path.basename(slc)
        dd = dd.replace(".vrt","") 
        dateList.append(dd) 

    dateList.sort()  
    return dateList

def getStackDict(dateList, miniStackDir, datumDir, outDir, miniStackSize, subDir = "EVD", fileExtension=".slc", outputExtension = ".phase"):

    miniStackDir = os.path.abspath(miniStackDir)
    datumDir = os.path.abspath(datumDir)
    outDir = os.path.abspath(outDir)

    stackSize = len(dateList)
    miniStackCount = 0;
    indStart = 0
    
    stackDict = {}
 
    while (indStart < stackSize):
        ##Update mini stack counter
        miniStackCount = miniStackCount + 1

        ###Determine end of ministack
        indEnd = indStart + inps.miniStackSize
        indEnd = min(indEnd, stackSize)
            
        miniStackDates = dateList[indStart:indEnd]
        
        miniStackPath = dateList[indStart] + "_" + dateList[indEnd-1]
        miniStackPath = os.path.join(miniStackDir, miniStackPath)
        #datumPath = os.path.join(datumDir , "EVD/" + dateList[indEnd-1]+".slc")
        datumPath = os.path.join(datumDir , subDir , dateList[indEnd-1]+ fileExtension)
        print(miniStackPath)
        for dd in miniStackDates:
            #miniStackDate = os.path.join(miniStackPath, "EVD/" + dd + ".slc")
            miniStackDate = os.path.join(miniStackPath, subDir, dd + fileExtension)
            temporalCoh = os.path.join(miniStackPath, "EVD/tcorr.bin")
            #outName = os.path.join(outDir, dd + "/" + dd + ".phase")
            outName = os.path.join(os.path.abspath(outDir), dd + outputExtension)
            stackDict[dd] = [miniStackDate, datumPath, temporalCoh, outName]

        indStart += inps.miniStackSize

    return stackDict

def getSize(data):

    ds = gdal.Open(data, gdal.GA_ReadOnly)
    length = ds.RasterYSize
    width = ds.RasterXSize
    ds = None
    return length, width

def adjustMiniStackPhase(miniStackSlc, adjustSlc, output, length, width):

    outDir = os.path.dirname(output)
    if not os.path.exists(outDir):
        os.makedirs(outDir)

    miniSlc = np.memmap(miniStackSlc, dtype=np.complex64, mode='r', shape=(length,width))
    adjust = np.memmap(adjustSlc, dtype=np.complex64, mode='r', shape=(length,width))
    outPhase = np.memmap(output, dtype=np.float32, mode='w+', shape=(length,width))
  
    ind = [ii for ii in range(0, length, 1000)]
    if ind[-1] != length:
        ind.append(length)

    print("writing to >>> ", output)

    for ii in range(len(ind)-1):
        startLine = ind[ii]
        endLine = ind[ii+1]
        miniSlc_block = miniSlc[startLine:endLine, :]  
        adjust_block = adjust[startLine:endLine, :]
        outPhase[startLine:endLine, :] = np.angle(miniSlc_block*np.conjugate(adjust_block))
    
    miniSlc = None
    adjust = None
    outPhase = None

    write_xml(output, length, width)

    return None

def write_xml(outFile, length, width):

    outImage = isceobj.Image.createImage()
    outImage.setFilename(outFile)
    outImage.setWidth(width)
    outImage.setLength(length)
    #unwImage.imageType = 'unw'
    outImage.bands = 1
    outImage.scheme = 'BIL'
    outImage.dataType = 'FLOAT'
    outImage.setAccessMode('read')
   # unwImage.createImage()
    outImage.renderHdr()
    outImage.renderVRT()

def adjust_acquisition_vrt(miniStackSlc, adjustSlc, output, length, width):

    vrttmpl = '''
<VRTDataset rasterXSize="{width}" rasterYSize="{length}">
  <VRTRasterBand dataType="Float32" band="1" subClass="VRTDerivedRasterBand">
    <Description>Sum</Description>
    <PixelFunctionType>sum</PixelFunctionType>
    <SimpleSource>
      <SourceFilename>{miniStackSlc}</SourceFilename>
    </SimpleSource>
    <SimpleSource>
      <SourceFilename >{adjustment}</SourceFilename>
    </SimpleSource>
  </VRTRasterBand>
</VRTDataset>'''
    
    with open( '{0}.vrt'.format(output) , 'w') as fid:
            fid.write( vrttmpl.format(width=width,
                                     length=length,
                                     miniStackSlc=miniStackSlc,
                                     adjustment=adjustSlc))
    
def adjust_acquisition_wrapped_vrt(miniStackSlc, adjustSlc, output, length, width):

    vrttmpl = '''
<VRTDataset rasterXSize="{width}" rasterYSize="{length}">
  <VRTRasterBand dataType="CFloat32" band="1" subClass="VRTDerivedRasterBand">
    <Description>multiply</Description>
    <PixelFunctionType>mul</PixelFunctionType>
    <SimpleSource>
      <SourceFilename>{miniStackSlc}</SourceFilename>
    </SimpleSource>
    <SimpleSource>
      <SourceFilename>{adjustment}</SourceFilename>
    </SimpleSource>
  </VRTRasterBand>
</VRTDataset>'''

    with open( '{0}.vrt'.format(output) , 'w') as fid:
            fid.write( vrttmpl.format(width=width,
                                     length=length,
                                     miniStackSlc=miniStackSlc,
                                     adjustment=adjustSlc))


def get_ref_phase(unwrap_file, ref_y, ref_x):
    ds = gdal.Open(unwrap_file+".vrt", gdal.GA_ReadOnly)
    b = ds.GetRasterBand(1)
    ref_phase = b.ReadAsArray(ref_x, ref_y, 1, 1)
    ds = None
    return -1*ref_phase[0][0]


def adjust_wrapped(dateList, inps):

    stackDict = getStackDict(dateList, inps.miniStackDir, inps.datumDir, inps.outDir, inps.miniStackSize, subDir ="EVD", fileExtension=".slc", outputExtension=".slc")

    length, width = getSize(stackDict[dateList[0]][0])
    print("length, width: {0}, {1}".format(length, width))

    for k in stackDict.keys():
            print(k)
            d1 = datetime.datetime(*time.strptime(k,"%Y%m%d")[0:6])
            day_of_year = d1.timetuple().tm_yday
            acq = float(d1.year)+float(day_of_year-1)/365.25

            miniStackSlc = stackDict[k][0]
            adjustSlc = stackDict[k][1]
            output = stackDict[k][3]
            print("mini stack slc: ", miniStackSlc)
            print("datum compensation slc: ", adjustSlc)
            print("adjusted output phase: ", output)
            print("*****************")
            adjust_acquisition_wrapped_vrt(miniStackSlc, adjustSlc, output, length, width)

def adjust_unwrapped(dateList, inps):

    stackDict = getStackDict(dateList, inps.miniStackDir, inps.datumDir 
                             , inps.outDir, inps.miniStackSize, subDir = "unwrap"
                             , fileExtension=".unw", outputExtension = ".phase")

    timeSeriesTemplate = '''

<VRTRasterBand dataType="Float32">
        <SimpleSource>
            <SourceFilename>{path}</SourceFilename>
        </SimpleSource>
        <Metadata domain="slc">
            <MDI key="Date">{date}</MDI>
            <MDI key="AcquisitionTime">{acq}</MDI>
        </Metadata>
        <ScaleOffset>{reference_phase}</ScaleOffset>
</VRTRasterBand>\n'''

    with open(inps.timeSeriesName, 'w') as fid_timeseries:
        fid_timeseries.write( '<VRTDataset rasterXSize="{xsize}" rasterYSize="{ysize}">\n'.format(xsize=width, ysize=length))
        for k in stackDict.keys():
            print(k)
            d1 = datetime.datetime(*time.strptime(k,"%Y%m%d")[0:6])
            day_of_year = d1.timetuple().tm_yday
            acq = float(d1.year)+float(day_of_year-1)/365.25

            miniStackSlc = stackDict[k][0]
            adjustSlc = stackDict[k][1]
            output = stackDict[k][3]
            print("mini stack slc: ", miniStackSlc)
            print("datum compensation slc: ", adjustSlc)
            print("adjusted output phase: ", output)
            print("*****************")
            adjust_acquisition_vrt(miniStackSlc+".vrt", adjustSlc+".vrt", output, length, width)
            reference_phase = get_ref_phase(output, inps.ref_y, inps.ref_x)

            fid_timeseries.write(timeSeriesTemplate.format(path=output+".vrt", date=k, acq=acq, reference_phase=reference_phase))

        fid_timeseries.write('</VRTDataset>')

        #adjustMiniStackPhase(miniStackSlc, adjustSlc, output, length, width)
    

if __name__ == '__main__':
    '''
    Main driver.
    '''
    #*************************************************************#
    # read the input options and prepare some output directories
    inps = cmdLineParser()
    if not os.path.exists(inps.outDir):
        os.makedirs(inps.outDir)

    # A list of acquisition dates
    dateList = getDates(inps.slcDir)

    if inps.unwrapped:
        adjust_unwrapped(dateList, inps)

    else:
        adjust_wrapped(dateList, inps)

    
