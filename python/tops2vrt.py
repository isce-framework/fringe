#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import numpy as np 
import os
import glob
from osgeo import gdal


def cmdLineParse():
    '''
    Command line parser.
    '''
    import argparse

    parser = argparse.ArgumentParser(description='Tops SLC stack to VRT')
    parser.add_argument('-i', '--input', dest='indir', type=str,
                        required=True,
                        help='Merged directory from ISCE2 Sentinel stack processor')
    parser.add_argument('-s', '--stack', dest='stackdir', type=str,
                        default='stack',
                        help='Directory where the co-registered SLC stack VRT will be stored (default: "stack")')
    parser.add_argument('-g', '--geom', dest='geomdir', type=str,
                        default='geometry',
                        help='Directory where the geometry VRTs will be stored (default: "geometry")')
    parser.add_argument('-b', '--bbox', dest='bbox', nargs=4, type=int,
                        default=None, metavar=('Y0', 'Y1', 'X0', 'X1'),
                        help='bounding box in row col: minLine maxLine minPixel maxPixel')
    parser.add_argument('-B', '--geo_bbox', dest='geobbox', nargs=4, type=float,
                        default=None, metavar=('S', 'N', 'W', 'E'),
                        help='bounding box in lat lon: South North West East')

    inps = parser.parse_args()

    return inps


def radarGeometryTransformer(latfile, lonfile, epsg=4326):
    '''
    Create a coordinate transformer to convert map coordinates to radar image line/pixels.
    '''

    driver = gdal.GetDriverByName('VRT')
    inds = gdal.OpenShared(latfile, gdal.GA_ReadOnly)
    tempds = driver.Create('', inds.RasterXSize, inds.RasterYSize, 0)
    inds = None

    tempds.SetMetadata({'SRS': 'EPSG:{0}'.format(epsg),
                        'X_DATASET': lonfile,
                        'X_BAND': '1',
                        'Y_DATASET': latfile,
                        'Y_BAND': '1',
                        'PIXEL_OFFSET': '0',
                        'LINE_OFFSET': '0',
                        'PIXEL_STEP': '1',
                        'LINE_STEP': '1'}, 'GEOLOCATION')

    trans = gdal.Transformer(tempds, None, ['METHOD=GEOLOC_ARRAY'])

    return trans


def lonlat2pixeline(lonFile, latFile, lon, lat):
    trans = radarGeometryTransformer(latFile, lonFile)

    ###Checkour our location of interest
    success, location = trans.TransformPoint(1, lon, lat, 0.)
    if not success:
        print('Location outside the geolocation array range')

    return location


def getLinePixelBbox(geobbox, latFile, lonFile):
    south, north, west, east = geobbox

    se = lonlat2pixeline(lonFile, latFile, east, south)
    nw = lonlat2pixeline(lonFile, latFile, west, north)

    ymin = int(np.round(np.min([se[1], nw[1]])))
    ymax = int(np.round(np.max([se[1], nw[1]])))

    xmin = int(np.round(np.min([se[0], nw[0]])))
    xmax = int(np.round(np.max([se[0], nw[0]])))

    print("x min-max: ", xmin, xmax)
    print("y min-max: ", ymin, ymax)

    return ymin, ymax, xmin, xmax


if __name__ == '__main__':
    '''
    Main driver.
    '''

    ##Parse command line
    inps = cmdLineParse()

    ###Get ann list and slc list
    slclist = glob.glob(os.path.join(inps.indir, 'SLC', '*', '*.slc.full.vrt'))
    num_slc = len(slclist)

    print('number of SLCs discovered: ', num_slc)
    slclist.sort()

    # Extract first slc width and height
    ds = gdal.Open(slclist[0], gdal.GA_ReadOnly)
    width = ds.RasterXSize
    height = ds.RasterYSize
    ds = None

    ####Set up single stack file
    if os.path.exists(inps.stackdir):
        print('stack directory: {0} already exists'.format(inps.stackdir))
    else:
        print('creating stack directory: {0}'.format(inps.stackdir))
        os.makedirs(inps.stackdir)

    latFile = os.path.join(inps.indir, "geom_reference", "lat.rdr.full.vrt")
    lonFile = os.path.join(inps.indir, "geom_reference", "lon.rdr.full.vrt")

    # setting up a subset of the stack
    if inps.geobbox:
        # if the bounding box in geo-coordinate is given, this has priority
        print("finding bbox based on geo coordinates of {} ...".format(
            inps.geobbox))
        ymin, ymax, xmin, xmax = getLinePixelBbox(inps.geobbox, latFile,
                                                  lonFile)

    elif inps.bbox:
        # if bbox in geo not given then look for line-pixel bbox
        print("using the input bbox based on line and pixel subset")
        ymin, ymax, xmin, xmax = inps.bbox

    else:
        # if no bbox provided, the take the full size
        ymin, ymax, xmin, xmax = [0, height, 0, width]

    xsize = xmax - xmin
    ysize = ymax - ymin

    slcs_base_file = os.path.join(inps.stackdir, 'slcs_base.vrt')
    print('write vrt file for stack directory')

    # Wavelength for Sentinel-1
    wavelength = 0.05546576
    with open(slcs_base_file, 'w') as fid:
        fid.write(
            '<VRTDataset rasterXSize="{xsize}" rasterYSize="{ysize}">\n'.format(
                xsize=xsize, ysize=ysize))

        for ind, slc in enumerate(slclist):
            # Extract SLCs width and height
            ds = gdal.Open(slc, gdal.GA_ReadOnly)
            width = ds.RasterXSize
            height = ds.RasterYSize
            ds = None

            # Extract date from co-registered SLCs
            date = os.path.basename(os.path.dirname(slc))

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
                                 date=date, acq=date,
                                 wvl=wavelength, index=ind + 1,
                                 path=slc)
            fid.write(outstr)

        fid.write('</VRTDataset>')

    ####Set up latitude, longitude and height files

    if os.path.exists(inps.geomdir):
        print('directory {0} already exists.'.format(inps.geomdir))
    else:
        print('creating geometry directory: {0}'.format(inps.geomdir))
        os.makedirs(inps.geomdir)

    vrttmpl = '''<VRTDataset rasterXSize="{xsize}" rasterYSize="{ysize}">
    <VRTRasterBand dataType="Float64" band="1">
      <SimpleSource>
        <SourceFilename>{PATH}</SourceFilename>
        <SourceBand>1</SourceBand>
        <SourceProperties RasterXSize="{width}" RasterYSize="{height}" DataType="Float64"/>
        <SrcRect xOff="{xmin}" yOff="{ymin}" xSize="{xsize}" ySize="{ysize}"/>
        <DstRect xOff="0" yOff="0" xSize="{xsize}" ySize="{ysize}"/>
      </SimpleSource>
    </VRTRasterBand>
</VRTDataset>'''

    print('write vrt file for geometry dataset')
    layers = ['lat', 'lon', 'hgt']
    for ind, val in enumerate(layers):
        with open(os.path.join(inps.geomdir, val + '.vrt'), 'w') as fid:
            fid.write(vrttmpl.format(xsize=xsize, ysize=ysize,
                                     xmin=xmin, ymin=ymin,
                                     width=width,
                                     height=height,
                                     PATH=os.path.abspath(
                                         os.path.join(inps.indir,
                                                      'geom_reference',
                                                      val + '.rdr.full.vrt')),
                                     linewidth=width * 8))
