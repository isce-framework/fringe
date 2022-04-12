#! /usr/bin/env python3

# Author: Marin Govorcin
# Caltech JPL, Dec 05, 2021

import os
import argparse
import xml.etree.ElementTree as ET
from xml.dom import minidom
import isce3
from isce3.unwrap import snaphu
from isce3.unwrap import Phass
from isce3.io.gdal import Raster as gRaster # Snaphu module takes isce3.io.gdal.Raster inputs
from isce3.io import Raster # Phass module takes isce3.io.Raster inputs

class ParseKwargs(argparse.Action):
    def __call__(self, parser, namespace, values, option_string=None):
        setattr(namespace, self.dest, dict())
        for value in values:
            key, value = value.split('=')
            getattr(namespace, self.dest)[key] = value

def cmdLineParser():
    
    '''
    Command Line Parser
    Script to use Phass/Snaphu in the isce3 environment

    USAGE:
    Phass
    unwrap_isce3.py -i 20190610_20190716.int -c tcorr_ds_ps.bin -k corr_thresh=0.3, good_corr=0.7
    unwrap_isce3.py -i 20190610_20190716.int -c tcorr_ds_ps.bin -k corr_thresh=0.3, good_corr=0.7 min_region=200

    Snaphu:
    unwrap_isce3.py -i 20190610_20190716.int -c tcorr_ds_ps.bin -u snaphu -k nlooks=27 tiles=[4,4,500,500] nproc=8 min_region=200
    '''
    unw_options = ''' phass  ::  corr_thresh=0.3, good_corr=0.7, min_region=200, \n
    snaphu :: nlooks=1, tiles=[1, 1, 0, 0], min_region=200, nproc=1, single_tile_opt=False '''

    parser = argparse.ArgumentParser(description='Phass Unwrapper', formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument('-i', '--interferogram_file', type=str, dest='igram',
                        required=True, help='Input interferogram file in complex format, example: file.int')
    parser.add_argument('-c', '--correlation_file', type=str, dest='correlation',
                        required=True, help='Input correlation file, example: tcorr_ds_ps.bin')
    parser.add_argument('-m', '--mask', type=str, dest='mask', default=None,
                        help='Input mask byte file')
    parser.add_argument('-o', '--output_dir', type=str, dest='outputDir',
                        default='./', help='Output directory path')
    parser.add_argument('-u', '--unwrapper', type=str, dest='unwrapper',
                        default='phass', help='Unwrapper method: phass/snaphu')
    parser.add_argument('-k', '--kwargs', type=str, dest='kwargs', nargs='*',
                        help=unw_options, action=ParseKwargs)

    return parser.parse_args()

############################################################################################################
############################ UTILITS FOR UNWRAPPING OUTPUT METDADATA #######################################
############################################################################################################

GDAL_DATATYPE = {
    1 : "Byte",
    2 : "UInt16",
    3 : "Int16",
    4 : "UInt32",
    5 : "UInt64",
    6 : "Float32",
    7 : "Float64",
    8 : "CInt16",
    9 : "CInt32",
    10 : "CFloat32",
    11 : "CFloat64"
    }

GDAL2ISCE = {
    1  : 'BYTE',
    2  : 'BYTE',
    3  : 'SHORT',
    5  : 'INT',
    6  : 'FLOAT',
    7  : 'DOUBLE',
    10 : 'CFLOAT',
    11 : 'CFLOAT'
    }

def write_vrt(infile, width, length, dtype):
    ## infile :: str  path to inputfile
    infile_dir  = os.path.dirname(os.path.abspath(infile))
    
    vrttmpl='''<VRTDataset rasterXSize="{width}" rasterYSize="{length}">
    <VRTRasterBand dataType="{DTYPE}" band="1" subClass="VRTRawRasterBand">
        <SourceFilename relativeToVRT="1">{PATH}</SourceFilename>
        <ImageOffset>0</ImageOffset>
        <PixelOffset>4</PixelOffset>
        <LineOffset>{linewidth}</LineOffset>
        <ByteOrder>LSB</ByteOrder>
    </VRTRasterBand>
</VRTDataset>'''
        
    outfile = os.path.join(infile_dir, '{0}.vrt'.format(infile))

    with open(outfile, 'w') as fid:
        fid.write(vrttmpl.format(width = width,
                                 length = length,
                                 DTYPE = GDAL_DATATYPE[dtype],
                                 PATH = infile,
                                 linewidth = 4 * width))

def property_xml(root, name, value, doc):
    property = ET.SubElement(root, 'property', {'name':name})
    ET.SubElement(property, 'value').text = value
    ET.SubElement(property, 'doc').text=doc

def component_xml(root, name, size, doc):
    coord = ET.SubElement(root, 'component', {'name':name})
    ET.SubElement(coord, 'factorymodule').text = 'isceobj.Image'
    ET.SubElement(coord, 'factoryname').text = 'createCoordinate'
    ET.SubElement(coord, 'doc').text = doc
    property_xml(coord, 'delta', '1.0', 'Coordinate quantization.')
    property_xml(coord, 'endingvalues', str(float(size)), 'Ending value of the coordinate.')
    property_xml(coord, 'family', 'imagecoordinate', 'Instance family name.')
    property_xml(coord, 'name', 'imagecoordinate_name', 'Instance name.')
    property_xml(coord, 'size', str(size), 'Coordinate size.')
    property_xml(coord, 'startingvalue', '0.0', 'Starting value of the coordinate.')

def write_xml(infile, width, length, bands, datatype, scheme):
    root = ET.Element('imageFile')
    property_xml(root, 'access_mode', 'READ', 'Image access mode.')
    property_xml(root, 'byte_order', 'l', 'Endianness of the image.')
    component_xml(root, 'coordinate1', width, 'First coordinate of a 2D image (width).')
    component_xml(root, 'coordinate2', length, 'Second coordinate of a 2D image (length).')
    property_xml(root, 'data_type', GDAL2ISCE[datatype], 'Image data type.')
    property_xml(root, 'family', 'image', 'Instance family name.')
    property_xml(root, 'file_name', os.path.abspath(infile), 'Name of the image file.')
    property_xml(root, 'length', str(length), 'Image length.')
    property_xml(root, 'name', 'image_name', 'Instance name.')
    property_xml(root, 'number_bands', str(bands), 'Number of image bands.')
    property_xml(root, 'scheme', scheme, 'Interleaving scheme of the image.')
    property_xml(root, 'width', str(width), 'Image width.')
    property_xml(root, 'x_max', str(width), 'Maximum range value.')
    property_xml(root, 'x_min', '0.0', 'Minimum range value.')
    
    #### Add ident to xml file
    xmlstr = minidom.parseString(ET.tostring(root)).toprettyxml(indent="   ")
    with open(infile+'.xml', "w") as f:
        f.write(xmlstr)

def unw_vrt_xml(infile, scheme):
    output = Raster(infile)
    width = output.width
    length = output.length
    dtype = output.datatype
    bands = output.band

    write_vrt(infile, width, length, dtype)
    write_xml(infile, width, length, bands, dtype, scheme)


############################################################################################################
###################################    UNWRAPPING PHASS/SNAPHU     #########################################
############################################################################################################

def vrt_cpx2rad(infile, width, length, stype):
    '''
    Function to extract interferogram phase from complex file
    needed when unwrapping interferogram in radians.
    
    Find the input interferogram dtype, if complex32 or complex64
    reate vrt file that extract phase values onthefly

    Input:
         infile  :: str  wrapped interferogram filename
    Output
         outfile :: str  output .vrt file
    '''
     
    # Get path to directory and name of inputfile
    infile_dir  = os.path.dirname(os.path.abspath(infile))
    infile_name = os.path.splitext(os.path.basename(infile))[0]

    #Input file is COMPLEX, create VRT with default pixelfunction : phase to extract phase in radian
    if stype == 10 or stype == 11:
        
        vrttmpl='''<VRTDataset rasterXSize="{width}" rasterYSize="{length}">
    <VRTRasterBand dataType="{DTYPE}" band="1" subClass="VRTDerivedRasterBand">
        <Description>Wrapped interferogram in radians on-the-fly</Description>
        <SimpleSource>
           <SourceFilename relativeToVRT="1">{PATH}</SourceFilename>
           <ImageOffset>0</ImageOffset>
           <PixelOffset>4</PixelOffset>
           <LineOffset>{linewidth}</LineOffset>
           <ByteOrder>LSB</ByteOrder>
        </SimpleSource>
        <PixelFunctionType>phase</PixelFunctionType>
        <SourceTransferType>{STYPE}</SourceTransferType>
    </VRTRasterBand>
</VRTDataset>'''
        
        if stype == 10:
            dtype = 6
        elif stype == 11:
            dtype = 7

        outfile = os.path.join(infile_dir, '{0}.phase.vrt'.format(infile_name))

        with open(outfile, 'w') as fid:
            fid.write(vrttmpl.format(width = width,
                                     length = length,
                                     DTYPE = GDAL_DATATYPE[dtype],
                                     STYPE = GDAL_DATATYPE[stype],
                                     PATH = infile,
                                     linewidth = 4 * width))

    else:
        outfile = os.path.abspath(infile)

    # Return the path to VRT file as input to Phass
    return outfile

def unwrap_phass(infile, corr_file, unw_filename, conncomp_filename,
                 correlation_threshold=0.3, good_correlation=0.7, min_region=200, print_msg=True):
    '''
    Unwrap wrapped interferogram [in radians] using Phass
    
    Input:
        infile                :: str    Wrapped phase in radians [interferogram *.int in PS_DS folder for fringe], input format .vrt
        corr_file             :: str    Correlation file [tcorr_ds_ps.bin for fringe], input format .vrt
        unw_filename          :: str    Output filename of unwrapped interferogram
        conncomp_filename     :: str    Output filename of connected components
        correlation_threshold :: float  Threshold for minimum coherence
        good_correlation      :: float  Good coherence
        min_region            :: int    Minimum number of pixels per region (connected component)

    Output:
        unw_igram   :: Raster, Unwrapped interferogram, type: isce3.io.Raster
        conn_comp   :: Raster, Connected component labels, type: isce3.io.Raster
    '''
    
    phass = Phass()
    
    # If input interferogram is in COMPLEX format, create .vrt file
    # to get interferogram phase in radians on-the-fly
    infile_name = infile
    infile = Raster(infile)
    phass_input = vrt_cpx2rad(infile_name, infile.width, infile.length, infile.datatype())

    # Input files
    igram = Raster(phass_input)
    correlation = Raster(corr_file)

    # Output files
    unw_igram = Raster(unw_filename, igram.width, igram.length,1, 6, "ENVI") #works also with ISCE format
    conn_comp = Raster(conncomp_filename, igram.width, igram.length, 1, 4, "ENVI")

    # Phass unwrapping parameters
    phass.correlation_threshold = correlation_threshold
    phass.good_correlation = good_correlation
    phass.min_pixels_region = min_region

    if print_msg:
        print('\nPhass Unwrapping:')
        print('Input:')
        print('  Phase: {}'.format(os.path.abspath(phass_input)))
        print('        width: {}, length {}, dtype: {}'.format(igram.width, igram.length, GDAL_DATATYPE[igram.datatype()]))
        print('  Correlation: {}'.format(os.path.abspath(corr_file)))
        print('        width: {}, length {}, dtype: {}'.format(correlation.width, correlation.length, GDAL_DATATYPE[correlation.datatype()]))
        print('  Parameters: correlation_threshold={}, good_correlation={}, min_region={}'.format(correlation_threshold, good_correlation, min_region))

    # Unwrap
    phass.unwrap(igram, correlation, unw_igram, conn_comp)

    if print_msg:
        print('#############################################')
        print('Finished unwrapping, output files:')
        print('  Unwrapped interferogram: {}'.format(os.path.abspath(unw_filename)))
        print('        width: {}, length {}, dtype: {}'.format(unw_igram.width, unw_igram.length, GDAL_DATATYPE[unw_igram.datatype()]))
        print('  Connected components: {}'.format(os.path.abspath(conncomp_filename)))
        print('        width: {}, length {}, dtype: {}'.format(conn_comp.width, conn_comp.length, GDAL_DATATYPE[conn_comp.datatype()]))

def unwrap_snaphu(infile, corr_file, unw_filename, conncomp_filename, mask=None, cost_mode='defo',
                  nlooks=1, tile_nrows=1, tile_ncols=1, row_overlap=0, col_overlap=0, tile_cost_thresh=500,
                  min_region=100,nproc=1, single_tile_opt=False, defo_thresh=1.25, print_msg=True):

    '''
    Unwrap wrapped interferogram using Snaphu
    
    Input:
        igram                 :: str    Input interferogram. Must be GDT_CFloat32 datatype
        corr_file             :: str    Correlation magnitude, normalized to the interval [0, 1], GDT_Float32 datatype
        unw_filename          :: str    Output filename of unwrapped interferogram
        conncomp_filename     :: str    Output filename of connected components
        mask                  :: str    Binary mask of valid pixels. Zeros in this raster indicate interferogram
                                        pixels that should be masked out, GDT_Byte datatype
        nlooks                :: float  Effective number of looks used to form the input correlation data
        tile_nrows            :: int    Number of tiles along the row directions
        tile_ncols            :: int    Number of tiles along the column directions
        row_overlap           :: int    Overlap in number of rows between neighboring tiles
        col_overlap           :: int    Overlap in number of columns between neighboring tiles
        tile_cost_thresh      :: int    Cost threshold to use for determining boundaries of reliable regions.
                                        Larger cost threshold implies smaller regions (safer, but more expensive computationally
        min_region            :: int    Minimum number of pixels per region (connected component)
        nproc                 :: int    Number of threads for parallel processing (use with tiles)
        single_tile_optv      :: bool   Single Tile Reoptimization
        defo_thresh           :: float  Factor applied to rho0 to get threshold for whether or not phase discontinuity is possible

    Output:
        unw_igram   :: Raster, Unwrapped interferogram, type: isce3.io.Raster
        conn_comp   :: Raster, Connected component labels, type: isce3.io.Raster
    '''

    # Input file
    igram = gRaster(infile)
    correlation = gRaster(corr_file)
    
    if mask:
        print('Using mask file {}'.format(mask))
        mask = gRaster(mask)

    # Output file
    unw_igram = gRaster(unw_filename, igram.width, igram.length, isce3.io.gdal.GDT_Float32, "ENVI") #works also with ISCE format
    conn_comp = gRaster(conncomp_filename, igram.width, igram.length, isce3.io.gdal.GDT_UInt32, "ENVI")
    
    # Snaphu parameters
    tiling_params = snaphu.TilingParams(nproc=nproc, tile_nrows=tile_nrows, tile_ncols=tile_ncols,
                            row_overlap=row_overlap, col_overlap=col_overlap, tile_cost_thresh=tile_cost_thresh,
                            min_region_size=min_region, single_tile_reoptimize=single_tile_opt)

    conncomp_params = snaphu.ConnCompParams()
    corr_bias_model_params = snaphu.CorrBiasModelParams(min_corr_factor = defo_thresh)

    cost_params = None
    solver_params = None
    #init_method = 'mcf' #not yet included in the isce3.unwrap.snaphu
    pwr = None
    unwest = None
    
    phase_stddev_model_params = None
    scratchdir = './'
    verbose = print_msg
    debug = False

    if print_msg:
        print('\nSnaphu Unwrapping:')
        print('Input:')
        print('  Phase: {}'.format(os.path.abspath(infile)))
        print('        width: {}, length {}, dtype: {}'.format(igram.width, igram.length, GDAL_DATATYPE[igram.datatype.value]))
        print('  Correlation: {}'.format(os.path.abspath(corr_file)))
        print('        width: {}, length {}, dtype: {}'.format(correlation.width, correlation.length, GDAL_DATATYPE[correlation.datatype.value]))
        print('  Parameters: nlooks={}, tiles={}, tile_cost_thresh={}, min_region={}, nproc={}, single_tile_opt={}, defo_thresh={}'.format(nlooks, tiles, tile_cost_thresh, min_region, nproc, single_tile_opt, defo_thresh))

    # Unwrap
    snaphu.unwrap(unw_igram, conn_comp, igram, correlation, nlooks,
                  "defo", cost_params, pwr, mask, unwest, tiling_params, solver_params,conncomp_params,
                   corr_bias_model_params, phase_stddev_model_params, scratchdir, verbose, debug)

    if print_msg:
        print('#############################################')
        print('Finished unwrapping, output files:')
        print('  Unwrapped interferogram: {}'.format(os.path.abspath(unw_filename)))
        print('        width: {}, length {}, dtype: {}'.format(unw_igram.width, unw_igram.length, GDAL_DATATYPE[unw_igram.datatype.value]))
        print('  Connected components: {}'.format(os.path.abspath(conncomp_filename)))
        print('        width: {}, length {}, dtype: {}'.format(conn_comp.width, conn_comp.length, GDAL_DATATYPE[conn_comp.datatype.value]))


####################################################################################
####################################################################################
if __name__ == '__main__':

    # Main driver

    inps = cmdLineParser()

    ######################### Prepare input and output files #######################
    input_name = os.path.splitext(os.path.basename(inps.igram))[0]
    unw_output = os.path.join(inps.outputDir, input_name + '.unw')
    conncomp_output = os.path.join(inps.outputDir, input_name + '_int.conn')

    # Need gdal vrt file to read correctly the tcorr_ds_ps.bin format
    corr = gRaster(inps.correlation)
    write_vrt(inps.correlation, corr.width, corr.length, corr.datatype.value)
    inps.correlation = os.path.abspath(inps.correlation + '.vrt')

    ################################ Set-up and Unwrap #############################
    if inps.unwrapper == 'phass':
        corr_thresh = 0.3
        good_corr = 0.7
        min_region = 200
        
        if inps.kwargs:
            if 'corr_thresh' in inps.kwargs:
                corr_thresh = float(inps.kwargs['corr_thresh'])
            if 'good_corr' in inps.kwargs:
                good_corr = float(inps.kwargs['good_corr'])
            if 'min_region' in inps.kwargs:
                min_region = int(inps.kwargs['min_region'])
           
        unwrap_phass(inps.igram + '.vrt', inps.correlation, unw_output, conncomp_output,
                     correlation_threshold=corr_thresh, good_correlation=good_corr, min_region=min_region)

    elif inps.unwrapper == 'snaphu':
        nlooks = 27
        tiles = [1,1,0,0]
        min_region = 200
        nproc = 1
        single_tile_opt = False
        defo_thresh = 1.25
        tile_cost_thresh = 500

        if inps.kwargs:
            if 'nlooks' in inps.kwargs:
                nlooks = int(inps.kwargs['nlooks'])
            if 'tiles' in inps.kwargs:
                tile = inps.kwargs['tiles'].strip('][').split(',')
                tiles = [int(i) for i in tile]
            if 'min_region' in inps.kwargs:
                min_region = int(inps.kwargs['min_region'])
            if 'nproc' in inps.kwargs:
                nproc = int(inps.kwargs['nproc'])
            if 'single_tile_opt' in inps.kwargs:
                single_tile_opt = int(inps.kwargs['single_tile_opt'])
            if 'defo_thresh' in inps.kwargs:
                defo_thresh = float(inps.kwargs['defo_thresh'])
            if 'tile_cost_thresh' in inps.kwargs:
                tile_cost_thresh = float(inps.kwargs['tile_cost_thresh'])

        unwrap_snaphu(inps.igram + '.vrt', inps.correlation, unw_output, conncomp_output, mask = inps.mask, nlooks=nlooks, 
                      tile_nrows=tiles[0], tile_ncols=tiles[1], row_overlap=tiles[2], col_overlap=tiles[3], min_region=min_region, 
                      nproc=nproc, single_tile_opt=single_tile_opt, defo_thresh=defo_thresh, tile_cost_thresh=tile_cost_thresh)

    ################################ Write output vrt and xml file #############################

    unw_vrt_xml(unw_output, "BIL")
    unw_vrt_xml(conncomp_output, "BIP")
