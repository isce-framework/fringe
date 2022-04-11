import os
from osgeo import gdal


def extractComplexSource(infile):
    '''
    Extract complex source tag from xml file.
    '''

    from xml.etree import ElementTree as ET
    with open(infile, 'r') as fid:
        instr = fid.read()

    root = ET.fromstring(instr)
    band = root.find('VRTRasterBand')
    src = band.find('ComplexSource')
    return ET.tostring(src, encoding='unicode')

def updateComplexSource(invrt, infile, xsize):
    '''
    Update source and band information.
    '''

    from xml.etree import ElementTree as ET

    with open(invrt, 'r') as fid:
        instr = fid.read()

    root = ET.fromstring(instr)
    bnd = root.find('VRTRasterBand')

    band = bnd.find('ComplexSource')
    ###Update filename
    src = band.find('SourceFilename')
    src.text = os.path.abspath(infile)

    ###Update band number
    num = band.find('SourceBand')
    num.text = "2"

    ###Update blocksize
    blk = band.find('SourceProperties')
    blk.attrib['BlockXSize'] = str(xsize)
    blk.attrib['BlockYSize'] = "1"

    outstr = ET.tostring(root, encoding='unicode')
    print(outstr)
    with open(invrt, 'w') as fid:
        fid.write(outstr)


class StackVRT(object):
    ''' A class to represent a stack using gdal VRT 
    '''
    def __init__(self, xSize=None, ySize=None, outname="stack.vrt"):
        self.stackName = outname

    def configure(self, sampleFile, band=1):
        ds = gdal.Open(sampleFile)
        self.xSize = ds.RasterXSize
        self.ySize = ds.RasterYSize
        self.DataType = ds.GetRasterBand(band).DataType
        ds = None

        self.driver = gdal.GetDriverByName('VRT')
        self.vrt = self.driver.Create(self.stackName, self.xSize, self.ySize, 0)  

        self.band = 0

    def _get_ref_phase(self, rasterFileName, ref_x, ref_y):
        ds = gdal.Open(rasterFileName, gdal.GA_ReadOnly)
        b = ds.GetRasterBand(1)
        ref_phase = b.ReadAsArray(ref_x, ref_y, 1, 1)
        ds = None
        return -1*ref_phase[0][0]
 
    def _extract_band(self, rasterFileName, band):
        sourceDir = os.path.join(os.path.dirname(self.stackName), "source_bands")
        if not os.path.exists(sourceDir):
            os.makedirs(sourceDir)

        bndVRT = os.path.join(sourceDir, os.path.basename(rasterFileName))
        gdalTranslateOpts = gdal.TranslateOptions(format='VRT',
                                        bandList = [band])
        gdal.Translate(bndVRT, rasterFileName,  options=gdalTranslateOpts)
        return os.path.abspath(bndVRT)

    def _setReference(self, rasterFile, originalRasterFile, reference):

        if reference is not None:
            ref_x = reference[0]
            ref_y = reference[1]
            ref_phase = self._get_ref_phase(rasterFile, ref_x, ref_y)
        else:
            ref_phase = 0.0

        ds = gdal.Open(rasterFile, gdal.GA_Update)
        ds.GetRasterBand(1).SetOffset(ref_phase)
        ds = None

        ##Translate to complex source
        
        tempVRT2 = rasterFile+".temp"
        opts = gdal.TranslateOptions(format='VRT', bandList=[1], unscale=True)
        gdal.Translate(tempVRT2, rasterFile, options=opts)
        updateComplexSource(tempVRT2, originalRasterFile, self.xSize)
        return tempVRT2

    def addBand(self, rasterFileName, sourceband = 1, metadata=None, reference = None):
        self.band+=1
        tempVRT = self._extract_band(rasterFileName, sourceband)
        tempVRT2 = self._setReference(tempVRT, rasterFileName, reference)
        self.vrt.AddBand(self.DataType)
        band = self.vrt.GetRasterBand(self.band)
        complexstr = extractComplexSource(tempVRT2)
        band.SetMetadata({'source_0' : complexstr}, 'vrt_sources')
        if metadata is not None:
            for k in metadata.keys():
                band.SetMetadata({k:metadata[k]})

    def addBand_old(self, rasterFileName, sourceband = 1, metadata=None, reference = None):

        if sourceband>1:
           rasterFileName = self._extract_band(rasterFileName, sourceband)

        self.band+=1
        print(self.band, reference)
        if reference is not None:
            ref_x = reference[0]
            ref_y = reference[1]
            ref_phase = self._get_ref_phase(rasterFileName, ref_x, ref_y)
        else:
            ref_phase = 0.0



        options = ['subClass=VRTRasterBand',
                   'band={0}'.format(self.band),
                   'SourceFilename={0}'.format(rasterFileName)]

        self.vrt.AddBand(self.DataType, options)

        ds = gdal.Open(rasterFileName, gdal.GA_ReadOnly)
        
        self.vrt.AddBand(self.DataType, ds.GetRasterBand(1))
        self.vrt.GetRasterBand(self.band).SetOffset(ref_phase)
        if metadata is not None:
            for k in metadata.keys():
                self.vrt.GetRasterBand(self.band).SetMetadata({k:metadata[k]})

        print("band : ", self.band)

    def close(self):
        #gdalTranslateOpts = gdal.TranslateOptions(format='VRT', 
        #                                      width=outXSize, height=outYSize,
        #                                      srcWin=[0,0,outXSize*xlooks, outYSize*ylooks],
        #                                      noData=noData, resampleAlg=method)
        gdalTranslateOpts = gdal.TranslateOptions(format='VRT',
                                                  unscale=True)

       
        gdal.Translate(self.stackName, self.vrt,  options=gdalTranslateOpts)

   #def crop(self, x0, y0, width, length)       
   #def multiLook()
   

