#!/usr/bin/env python3

import glob
import StackVRT
import os

def cmdLineParse():
    '''
    Command line parser.
    '''
    import argparse

    parser = argparse.ArgumentParser(description='create a stack dataset using VRT')
    parser.add_argument('-i', '--input', dest='inputDir', type=str,
            required=True, help='The input parent directory (e.g., /u/k-data/dbekaert/APS_raytracing/Mexico/ALOS/track_188/processing/Igrams)')
    parser.add_argument('-p', '--pattern', dest='pattern', type=str,
            required=True, help='A name pattern which when combined with the input directory, lists all input files. e.g., "*/filt*conncomp"')
    parser.add_argument('-s', '--stack_name', dest='stackname', type=str,
            default='stack.vrt', help='full path of the output stack vrt file. (default is "./stack.vrt")')
    parser.add_argument('-r', '--reference_pixel', dest='reference_pixel', type=int, nargs="+",
            default=None, help='reference pixel [x,y]; e.g. -r 240 2145')

    parser.add_argument('-b', '--bbox', dest='bbox', nargs='+' , type=int, default=None,
            help='bounding box : minLine maxLine minPixel maxPixel')

    inps = parser.parse_args()

    return inps

if __name__ == '__main__':
    '''
    Main driver.
    '''

    ##Parse command line
    inps = cmdLineParse()
    outDir = os.path.dirname(inps.stackname)
    if outDir:
       if not os.path.exists(outDir):
          os.makedirs(outDir)

    files = glob.glob(os.path.join(inps.inputDir, inps.pattern))
    print("Total number of files:", len(files))
    st = StackVRT.StackVRT(outname=inps.stackname)
    st.configure(files[0])
    for ff in files:
        reference_follower = os.path.basename(os.path.dirname(ff)) 
        metaData = {"Dates":reference_follower}
        st.addBand(ff, metadata=metaData, sourceband=2, reference=inps.reference_pixel)

    st.close()
    st = None


