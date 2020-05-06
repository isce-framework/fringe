#!/usr/bin/env python3

# Author: Heresh Fattahi

import os
import glob
import argparse
import numpy as np
import gdal
import isce
import isceobj

#from Stack import Stack, MiniStack

def cmdLineParser():
    '''
    Command line parser.
    '''

    parser = argparse.ArgumentParser(description = 'creates a run file for unwrapping ministack and datum adjutment phase series',
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
            required=False, help='output directory')

    parser.add_argument('-u', '--unwrapper_script', type=str, dest='unwrapperScript',
            default="unwrap_fringe.py", help='unwraper script. e.g., unwrap_fringe.py')

    parser.add_argument('--unw_method','--unwrap_method', type=str, dest='unwrapMethod',
            default='phass', choices=('snaphu','phass'), help='unwrap method.')

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

def getStackDict(dateList, miniStackDir, datumDir, outDir, miniStackSize):

    stackSize = len(dateList)
    miniStackCount = 0;
    indStart = 0
    
    stackDict = {}
    datumDict = {}

    while (indStart < stackSize):
        ##Update mini stack counter
        miniStackCount = miniStackCount + 1

        ###Determine end of ministack
        indEnd = indStart + inps.miniStackSize
        indEnd = min(indEnd, stackSize)
            
        miniStackDates = dateList[indStart:indEnd]
        
        miniStackPath = dateList[indStart] + "_" + dateList[indEnd-1]
        miniStackPath = os.path.join(miniStackDir, miniStackPath)
        datumPath = os.path.join(datumDir , "EVD/" + dateList[indEnd-1]+".slc")
        datumUnwrap = os.path.join(datumDir , "unwrap", dateList[indEnd-1]+".unw")
        datumCor = os.path.join(datumDir , "EVD/" + "tcorr.bin")
        datumDict[dateList[indEnd-1]] = [datumPath, datumCor, datumUnwrap]
        print(miniStackPath)
        for dd in miniStackDates:
            miniStackDate = os.path.join(miniStackPath, "EVD/" + dd + ".slc")
            temporalCoh = os.path.join(miniStackPath, "EVD/tcorr.bin")
            #outName = os.path.join(outDir, dd + "/" + dd + ".unw")
            outName = os.path.join(miniStackPath, "unwrap" , dd + ".unw")
            stackDict[dd] = [miniStackDate, temporalCoh, outName]
            #stackDict[dd] = [miniStackDate, datumPath, temporalCoh, outName]
        print("***********")
        #################################

        indStart += inps.miniStackSize

    return stackDict, datumDict

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

if __name__ == '__main__':
    '''
    Main driver.
    '''
    #*************************************************************#
    # read the input options and prepare some output directories
    inps = cmdLineParser()

    # A list of acquisition dates
    dateList = getDates(inps.slcDir)
    stackDict, datumDict = getStackDict(dateList, inps.miniStackDir, inps.datumDir, inps.outDir, inps.miniStackSize)

    length, width = getSize(stackDict[dateList[0]][0]) 
    print("length, width: {0}, {1}".format(length, width)) 

    run_outname = "run_unwrap.sh"
    runf= open(run_outname,'w') 

    for k in stackDict.keys():
        print(k)
        miniStackSlc = stackDict[k][0]
        coh = stackDict[k][1]
        output = stackDict[k][2]
        print("mini stack slc: ", miniStackSlc)
        print("coherence: ", coh) 
        print("adjusted output phase: ", output)
        print("*****************")
        cmd = inps.unwrapperScript + " -m " + inps.unwrapMethod + " -i " + miniStackSlc + " -c " + coh + " -o " + output
        runf.write(cmd + "\n") 

    for k in datumDict.keys():
        print(k)
        datumSlc = datumDict[k][0]
        coh = datumDict[k][1]
        output = datumDict[k][2]
        print("mini stack slc: ", datumSlc)
        print("coherence: ", coh)
        print("adjusted output phase: ", output)
        print("*****************")
        cmd = inps.unwrapperScript + " -m " + inps.unwrapMethod + " -i " + datumSlc + " -c " + coh + " -o " + output
        runf.write(cmd + "\n")

    runf.close()
 


