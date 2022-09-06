#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import numpy as np 
import os
import glob
import datetime
from osgeo import gdal

def cmdLineParse():
    '''
    Command line parser.
    '''
    import argparse

    parser = argparse.ArgumentParser(description='Stripmap SLC stack to VRT')
    parser.add_argument('-i', '--input', dest='indir', type=str,
            required=True, help='Directory with segment of UAVSAR data')
    parser.add_argument('-s', '--stack', dest='stackdir', type=str,
            default='stack', help='Directory with stack vrts')
    parser.add_argument('-g', '--geom', dest='geomdir', type=str,
            default='geometry', help='Directory with geometry vrts')
    parser.add_argument('-c', '--slcs', dest='outdir', type=str,
            default='slcs', help='Directory with individual slc vrts')

    parser.add_argument('-b', '--bbox', dest='bbox', nargs=4, type=int, default=None, metavar=('Y0','Y1','X0','X1'),
            help='bounding box : minLine maxLine minPixel maxPixel')

    inps = parser.parse_args()

    return inps

if __name__ == '__main__':
    '''
    Main driver.
    '''
    
    ##Parse command line
    inps = cmdLineParse()

    ###Get ann list and slc list
    slclist = glob.glob(os.path.join(inps.indir, '*/*.slc'))

    print('Number of SLCs discovered: ', len(slclist))
    print('We assume that the SLCs and the ann files are sorted in the same order')
    
    slclist.sort()


    ###Read the first ann file to get some basic things like dimensions
    ###Just walk through each of them and create a separate VRT first
    if not os.path.exists(inps.outdir):
        print('Creating directory: {0}'.format(inps.outdir))
        os.makedirs(inps.outdir)
    else:
        print('Directory {0} already exists.'.format(inps.outdir))

    data = []
    dates = []

    width = None
    height = None

    for ind, slc in enumerate(slclist):

        ###Parse the ann file.
        metadata = {}
        width = None
        height = None
        path = None

        ds = gdal.Open(slc + '.vrt', gdal.GA_ReadOnly)
        width = ds.RasterXSize
        height = ds.RasterYSize
        ds = None


        metadata['WAVELENGTH'] = 0.03 
        metadata['ACQUISITION_TIME'] = os.path.basename(os.path.dirname(slc))
        
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

#        outname =  datetime.datetime.strptime(tag.upper(), '%d-%b-%Y %H:%M:%S UTC').strftime('%Y%m%d')

        outname = metadata['ACQUISITION_TIME']
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

    # setting up a subset of the stack
    ymin, ymax, xmin, xmax = [0 , height, 0 , width]
    if inps.bbox:
        ymin, ymax, xmin, xmax = inps.bbox

    xsize = xmax - xmin
    ysize = ymax - ymin

    with open( os.path.join(inps.stackdir, 'slcs_base.vrt'), 'w') as fid:
        #fid.write( '<VRTDataset rasterXSize="{width}" rasterYSize="{height}">\n'.format(width=width, height=height))
        fid.write( '<VRTDataset rasterXSize="{xsize}" rasterYSize="{ysize}">\n'.format(xsize=xsize, ysize=ysize))

        for ind, (date, meta) in enumerate( zip(dates, data)):
            outstr = '''    <VRTRasterBand dataType="CFloat32" band="{index}">
        <SimpleSource>
            <SourceFilename>{path}</SourceFilename>
            <SourceBand>1</SourceBand>
            <SourceProperties RasterXSize="{width}" RasterYSize="{height}" DataType="CFloat32"/>
            <SrcRect xOff="{xmin}" yOff="{ymin}" xSize="{xsize}" ySize="{ysize}"/>
            <DstRect xOff="0" yOff="0" xSize="{xsize}" ySize="{ysize}"/>
        </SimpleSource>
        <Metadata domain="slc">
            <MDI key="Date">{date}</MDI>
            <MDI key="Wavelength">{wvl}</MDI>
            <MDI key="AcquisitionTime">{acq}</MDI>
        </Metadata>
    </VRTRasterBand>\n'''.format(width=width, height=height,
                                xmin=xmin, ymin=ymin,
                                xsize=xsize, ysize=ysize,
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


    
    vrttmpl='''<VRTDataset rasterXSize="{xsize}" rasterYSize="{ysize}">
    <VRTRasterBand dataType="Float64" band="1" subClass="VRTRawRasterBand">
        <SourceFilename>{PATH}</SourceFilename>
        <ImageOffset>0</ImageOffset>
        <PixelOffset>8</PixelOffset>
        <LineOffset>{linewidth}</LineOffset>
        <ByteOrder>LSB</ByteOrder>
        <SrcRect xOff="{xmin}" yOff="{ymin}" xSize="{xsize}" ySize="{ysize}"/>
        <DstRect xOff="0" yOff="0" xSize="{xsize}" ySize="{ysize}"/>
    </VRTRasterBand>
</VRTDataset>'''
    

    layers = ['lat', 'lon', 'z']

    for ind, val in enumerate(layers):
        with open( os.path.join(inps.geomdir, val+'.vrt'), 'w') as fid:
            fid.write( vrttmpl.format( xsize = xsize, ysize = ysize,
                                       xmin = xmin, ymin = ymin,
                                       width = width,
                                       height = height,
                                       PATH = os.path.abspath( os.path.join(inps.indir, 'geom_reference', val+'.rdr')),
                                       linewidth = width * 8))



