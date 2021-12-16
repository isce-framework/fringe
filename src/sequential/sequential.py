#!/usr/bin/env python3

# Author: Heresh Fattahi

import os
import argparse
import shutil
from Stack import Stack, MiniStack

def cmdLineParser():
    '''
    Command line parser.
    '''

    parser = argparse.ArgumentParser(description = 'Perform MLE-based phase-linking on a stack of coregistered SLCs',
                formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('-i', '--inDir', type=str, dest='inputDir',
            required=True, help='Input folder which contains folders for each SLC')
    parser.add_argument('-w', '--weight_dataset', type=str, dest='weightDS',
            required=True, help='Input weights dataset')
    parser.add_argument('-o', '--outDir', type=str, dest='outputDir',
            required=True, help='Output folder')
    parser.add_argument('-l', '--linesperblock', type=int, dest='linesPerBlock',
            default=64, help='Quantum for block of lines')
    parser.add_argument('-r', '--ram', type=int, dest='memorySize',
            default=2048, help='Memory in Mb to use')
    parser.add_argument('-x', '--xhalf', type=int, dest='halfWindowX',
            default=29, help='Half window size (range)')
    parser.add_argument('-y', '--yhalf', type=int, dest='halfWindowY',
            default=9, help='Half window size (azimuth)')
    parser.add_argument('-m', '--minneigh', type=int, dest='minNeighbors',
            default=5, help='Minimum number of neighbors for computation')

    parser.add_argument('-b', '--bbox', dest='bbox', nargs='+' , type=str, default=None,
            help='bounding box : minLine maxLine minPixel maxPixel \nORcoreg_stack/slcs_base.vrt')
    parser.add_argument('-s', '--mini_stack_size', type=int, dest='miniStackSize',
            default=10, help='mini stack size')

    parser.add_argument('-f', '--force', dest='forceprocessing', action='store_true',
            default=False, help='Force reprocessing')

    return parser.parse_args()


def vrt_file2bbox(vrt_file):
    '''
    Grab bounding box info from VRT file
    '''
    import defusedxml.ElementTree as ET
    print('read bbox info from VRT file: {}'.format(vrt_file))

    root = ET.parse(vrt_file).getroot()
    type_tag = root.find('VRTRasterBand/SimpleSource/SrcRect')
    xmin = int(type_tag.get('xOff'))
    ymin = int(type_tag.get('yOff'))
    xsize = int(type_tag.get('xSize'))
    ysize = int(type_tag.get('ySize'))
    xmax = xmin + xsize
    ymax = ymin + ysize

    bbox = (ymin, ymax, xmin, xmax)
    return bbox


def runEvd(inps, inputDataset, weightDS, outDir, miniStackCount, compressedSlcDir=None, compressedSlcName=None):
    '''
    Actually run Evd.
    '''

    import os
    os.environ['OPENBLAS_NUM_THREADS'] = "1"

    import evdlib
    aa = evdlib.Evd()
    
    ###Explicit wiring. Can be automated later.
    aa.inputDS = inputDataset
    aa.weightsDS = weightDS
    aa.outputFolder = outDir
    
    aa.miniStackCount = miniStackCount

    aa.blocksize = inps.linesPerBlock
    aa.memsize = inps.memorySize
    aa.halfWindowX = inps.halfWindowX
    aa.halfWindowY = inps.halfWindowY
    aa.minimumNeighbors = inps.minNeighbors

    if compressedSlcDir:
        aa.outputCompressedSlcFolder = compressedSlcDir
    else:
        aa.outputCompressedSlcFolder = aa.outputFolder

    if compressedSlcName:
        aa.compSlc = compressedSlcName
    else:
        aa.compSlc = "compslc.bin"

    aa.run()

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
    aa.run()


def moveCompressedSlc(outDir, compressedSlcDir):

    compSlc = os.path.join(outDir , "compslc.bin")
    compSlcHdr = os.path.join(outDir , "compslc.bin.hdr")

    tempCompSlcDir = os.path.join(compressedSlcDir , endDate)
    if not os.path.exists(tempCompSlcDir):
        os.makedirs(tempCompSlcDir)

    compSlcMoved = os.path.join(tempCompSlcDir, "compslc.slc") 
    compSlcHdrMoved = os.path.join(tempCompSlcDir, "compslc.slc.hdr")

    shutil.move(compSlc, compSlcMoved)
    shutil.move(compSlcHdr, compSlcHdrMoved)
    os.system("gdal_translate -of VRT {0} {1}".format(compSlcMoved , compSlcMoved + ".vrt"))

    return None

if __name__ == '__main__':
    '''
    Main driver.
    '''
    #*************************************************************#
    # read the input options and prepare some output directories
    inps = cmdLineParser()
    weightDS = inps.weightDS
    inps.outputDir = os.path.abspath(inps.outputDir)
    outDir = os.path.join(inps.outputDir, "fullStack")
    compressedSlcDir = os.path.join(inps.outputDir, "compressedSlc")
    if not os.path.exists(compressedSlcDir):
       os.makedirs(compressedSlcDir)

    #*************************************************************#
    # get vrt files for whole stack

    stack = Stack(inps.inputDir)

    # read bounding box input
    if inps.bbox is not None:
        if os.path.isfile(inps.bbox[0]):
            inps.bbox = vrt_file2bbox(inps.bbox[0])
        else:
            inps.bbox = tuple([int(i) for i in inps.bbox])
        print('input bounding box in (y0, y1, x0, x1): {}'.format(inps.bbox))
    stack.bbox = inps.bbox

    stack.gatherSLCs()
    stack.getDates()
    stack.configure(outDir)
    stack.writeStackVRT()

    # run calamp, ampdispersion, nmap to get wights (Weights won't change any )
    # This weight dataset will be used for all mini stacks
    #For testing
    #********************************#

    # get vrt files for first mini stack

    runEvdOnMiniStack = True


    
    #*************************************************************#
    # while there are more mini stacks kepp process them
    miniStackCount = 0;
    indStart = 0

    while (indStart < stack.size):
        ##Update mini stack counter
        miniStackCount = miniStackCount + 1

        ###Determine end of ministack
        indEnd = indStart + inps.miniStackSize
        indEnd = min(indEnd, stack.size)

        ###Figure out folder names
        startDate = stack.dateList[indStart]
        endDate = stack.dateList[indEnd-1]
        outDir = os.path.join(inps.outputDir, "miniStacks/" + startDate + "_" + endDate)

        if os.path.isdir(outDir) and (not inps.forceprocessing):
            print('{0} looks like it has already been processed. Skipping ... '.format(outDir))
        else:
            if (inps.forceprocessing):
                print('User has asked for reprocessing - {0}'.format(outDir))
            print('Processing {0}'.format(outDir))

            miniStack = MiniStack()
            miniStack.slcList = stack.slcList[indStart:indEnd]
            miniStack.getSize()
            miniStack.applyBbox = [True]*miniStack.size
            miniStack.updateMiniStack(compressedSlcDir)
            miniStack.bbox = inps.bbox
            miniStack.getDates()
            miniStack.configure(outDir)
            miniStack.writeStackVRT()

            outDir = os.path.join(outDir, "EVD")
            compressedSlcName = miniStack.dateList[-1] + ".slc"
            compSlcDir = os.path.join(compressedSlcDir, miniStack.dateList[-1])
            if not os.path.exists(compSlcDir):
                os.makedirs(compSlcDir)
            
            runEvd(inps, miniStack.stackVRT, weightDS, outDir, miniStackCount, compSlcDir, compressedSlcName)
            
            os.system("gdal_translate -of VRT {0} {1}".format(os.path.join(compSlcDir, compressedSlcName) , os.path.join(compSlcDir, compressedSlcName) + ".vrt"))

            # move compressed SLC of the current ministack to the compressed folder
            #moveCompressedSlc(outDir, compressedSlcDir)

            del miniStack


        ###Update index of ministack start
        indStart += inps.miniStackSize

     
    #*************************************************************#
    # Datum connection to connect the phase series of mini stacks
    # 
    outDir = os.path.join(inps.outputDir,"Datum_connection")     
    compSlcStack = Stack(compressedSlcDir)
    compSlcStack.gatherSLCs()
    compSlcStack.bbox = None
    compSlcStack.applyBbox = [False]*compSlcStack.size
    compSlcStack.getDates()
    compSlcStack.configure(outDir)
    compSlcStack.writeStackVRT()
    runEvd(inps, compSlcStack.stackVRT, weightDS, outDir + "/EVD", 1)
 
       

