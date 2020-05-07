#!/usr/bin/env python3
# -*- coding: utf-8 -*-

try:
    import isce
except ImportError:
    raise Exception('Cannot read in uavsar annotation files without ISCE. Please install ISCE before proceeding. If you have already installed ISCE, make sure that the PATHs and PYTHONPATHs are set correctly for using it.')

import numpy as np 
import os
import glob
from iscesys.Parsers import rdf
import datetime

def cmdLineParse():
    '''
    Command line parser.
    '''
    import argparse

    parser = argparse.ArgumentParser(description='UAVSAR SLC stack to VRT')
    parser.add_argument('-i', '--input', dest='indir', type=str,
            required=True, help='Directory with segment of UAVSAR data')
    parser.add_argument('-s', '--stack', dest='stackdir', type=str,
            default='stack', help='Directory with stack vrts')
    parser.add_argument('-g', '--geom', dest='geomdir', type=str,
            default='geometry', help='Directory with geometry vrts')
    parser.add_argument('-c', '--slcs', dest='outdir', type=str,
            default='slcs', help='Directory with individual slc vrts')

    inps = parser.parse_args()

    return inps

if __name__ == '__main__':
    '''
    Main driver.
    '''
    
    ##Parse command line
    inps = cmdLineParse()

    ###Get ann list and slc list
    slclist = glob.glob(os.path.join(inps.indir, '*1x1.slc'))
    annlist = glob.glob(os.path.join(inps.indir, '*.ann'))

    if len(slclist) != len(annlist):
        raise Exception('Number of SLCs {0} and Number of ann files {1} dont match'.format(len(slclist), len(annlist)))
    
    print('Number of SLCs discovered: ', len(slclist))
    print('We assume that the SLCs and the ann files are sorted in the same order')
    
    slclist.sort()
    annlist.sort()


    ###Read the first ann file to get some basic things like dimensions
    ###Just walk through each of them and create a separate VRT first

    if not os.path.exists(inps.outdir):
        print('Creating directory: {0}'.format(inps.outdir))
        os.makedirs(inps.outdir)
    else:
        print('Directory {0} already exists.'.format(inps.outdir))

    data = []
    dates = []

    for ind, (ann,slc) in enumerate( zip(annlist, slclist)):

        ###Parse the ann file.
        metadata = {}
        width = None
        height = None
        path = None

        parser = rdf.parse(ann)
    
        width = int(parser['slc_1_1x1 Columns'])
        height = int(parser['slc_1_1x1 Rows'])
        metadata['WAVELENGTH'] = parser['Center Wavelength']
        metadata['ACQUISITION_TIME'] = parser['Start Time of Acquisition']
        path = os.path.abspath(slc)

        tag = metadata['ACQUISITION_TIME'] 

        vrttmpl='''<VRTDataset rasterXSize="{width}" rasterYSize="{height}">
    <VRTRasterBand dataType="CFloat32" band="1" subClass="VRTRawRasterBand">
        <sourceFilename>{PATH}</sourceFilename>
        <ImageOffset>0</ImageOffset>
        <PixelOffset>8</PixelOffset>
        <LineOffset>{linewidth}</LineOffset>
        <ByteOrder>LSB</ByteOrder>
    </VRTRasterBand>
</VRTDataset>'''

        outname =  datetime.datetime.strptime(tag.upper(), '%d-%b-%Y %H:%M:%S UTC').strftime('%Y%m%d')
        with open( os.path.join(inps.outdir, '{0}.vrt'.format(outname)) , 'w') as fid:
            fid.write( vrttmpl.format(width=width,
                                     height=height,
                                     PATH=path,
                                     linewidth=8*width))

        
        data.append(metadata)
        dates.append(outname)


    ####Set up single stack file
    if os.path.exists( inps.stackdir):
        print('Stack directory: {0} already exists'.format(inps.stackdir))
    else:
        print('Creating stack directory: {0}'.format(inps.stackdir))
        os.makedirs(inps.stackdir)

    with open( os.path.join(inps.stackdir, 'slcs_base.vrt'), 'w') as fid:
        fid.write( '<VRTDataset rasterXSize="{width}" rasterYSize="{height}">\n'.format(width=width, height=height))

        for ind, (date, meta) in enumerate( zip(dates, data)):
            outstr = '''    <VRTRasterBand dataType="CFloat32" band="{index}">
        <SimpleSource>
            <SourceFilename>{path}</SourceFilename>
            <SourceBand>1</SourceBand>
            <SourceProperties RasterXSize="{width}" RasterYSize="{height}" DataType="CFloat32"/>
            <SrcRect xOff="0" yOff="0" xSize="{width}" ySize="{height}"/>
            <DstRect xOff="0" yOff="0" xSize="{width}" ySize="{height}"/>
        </SimpleSource>
        <Metadata domain="slc">
            <MDI key="Date">{date}</MDI>
            <MDI key="Wavelength">{wvl}</MDI>
            <MDI key="AcquisitionTime">{acq}</MDI>
        </Metadata>
    </VRTRasterBand>\n'''.format(width=width, height=height,
                                date=date, acq=meta['ACQUISITION_TIME'],
                                wvl = meta['WAVELENGTH'], index=ind+1, 
                                path = os.path.abspath( os.path.join(inps.outdir, date+'.vrt')))
            fid.write(outstr)

        fid.write('</VRTDataset>')

    ####Set up latitude, longitude and height files
    
    if os.path.exists( inps.geomdir):
        print('Directory {0} already exists.'.format(inps.geomdir))
    else:
        print('Creating geometry directory: {0}'.format(inps.geomdir))
        os.makedirs( inps.geomdir)


    llhfile = glob.glob(os.path.join(inps.indir, '*.llh'))
    if len(llhfile) == 0:
        raise Exception('No LLH file found')

    if len(llhfile) > 1:
        raise Exception('Multiple LLH files found')

    llhlookstag = os.path.splitext(llhfile[0])[0][-5:]
    
    mlwidth = parser['llh_{0} Columns'.format(llhlookstag)]
    mlheight = parser['llh_{0} Rows'.format(llhlookstag)]

    vrttmpl='''<VRTDataset rasterXSize="{width}" rasterYSize="{height}">
    <VRTRasterBand dataType="Float32" band="1" subClass="VRTRawRasterBand">
        <SourceFilename>{PATH}</SourceFilename>
        <ImageOffset>{ioff}</ImageOffset>
        <PixelOffset>12</PixelOffset>
        <LineOffset>{linewidth}</LineOffset>
        <ByteOrder>LSB</ByteOrder>
    </VRTRasterBand>
</VRTDataset>'''


    vrtstretchtmpl='''<VRTDataset rasterXSize="{width}" rasterYSize="{height}">
    <VRTRasterBand dataType="Float32" band="1">
        <SimpleSource resampling="bilinear">
            <SourceFilename>{src}</SourceFilename>
            <SourceBand>1</SourceBand>
            <SourceProperties RasterXSize="{mlwidth}" RasterYSize="{mlheight}" DataType="Float32"/>
            <SrcRect xOff="0" yOff="0" xSize="{mlwidth}" ySize="{mlheight}"/>
            <DstRect xOff="0" yOff="0" xSize="{width}" ySize="{height}"/>
        </SimpleSource>
    </VRTRasterBand>
</VRTDataset>'''



    layers = ['lat', 'lon', 'hgt']

    for ind, val in enumerate(layers):
        with open( os.path.join(inps.geomdir, '{1}_{0}.vrt'.format(llhlookstag, val)) , 'w') as fid:
            fid.write( vrttmpl.format(  width = mlwidth,
                                    height = mlheight,
                                    PATH = os.path.abspath(llhfile[0]),
                                    ioff=ind*4,
                                    linewidth = 12*mlwidth))

        with open( os.path.join(inps.geomdir, val+'.vrt'), 'w') as fid:
            fid.write( vrtstretchtmpl.format( width = width,
                                          height = height,
                                          mlwidth = mlwidth,
                                          mlheight = mlheight,
                                          src = os.path.abspath( os.path.join(inps.geomdir, '{1}_{0}.vrt'.format(llhlookstag, val) )) ))



    ###Setting up look vector
    lkvfile = glob.glob(os.path.join(inps.indir, '*.lkv'))
    if len(lkvfile) == 0:
        raise Exception('No LKV file found')

    if len(lkvfile) > 1:
        raise Exception('Multiple LKV files found')
    
    lkvlookstag = os.path.splitext(lkvfile[0])[0][-5:]
    

    mlwidth = parser['llh_{0} Columns'.format(llhlookstag)]
    mlheight = parser['llh_{0} Rows'.format(llhlookstag)]


    layers = ['elook', 'nlook', 'ulook']

    for ind, val in enumerate(layers):
        with open( os.path.join(inps.geomdir, '{1}_{0}.vrt'.format(lkvlookstag, val)) , 'w') as fid:
            fid.write( vrttmpl.format(  width = mlwidth,
                                    height = mlheight,
                                    PATH = os.path.abspath(lkvfile[0]),
                                    ioff=ind*4,
                                    linewidth = 12*mlwidth))

        with open( os.path.join(inps.geomdir, val+'.vrt'), 'w') as fid:
            fid.write( vrtstretchtmpl.format( width = width,
                                          height = height,
                                          mlwidth = mlwidth,
                                          mlheight = mlheight,
                                          src = os.path.abspath( os.path.join(inps.geomdir, '{1}_{0}.vrt'.format(lkvlookstag, val) )) ))

