#!/usr/bin/env python3

import os
import glob
from osgeo import gdal


vrttmpl='''
<VRTDataset rasterXSize="{width}" rasterYSize="{height}">
    <VRTRasterBand dataType="CFloat32" band="1" subClass="VRTRawRasterBand">
        <sourceFilename>{PATH}</sourceFilename>
        <ImageOffset>0</ImageOffset>
        <PixelOffset>8</PixelOffset>
        <LineOffset>{linewidth}</LineOffset>
        <ByteOrder>LSB</ByteOrder>
    </VRTRasterBand>
</VRTDataset>'''


vrttmpl2 = '''
    <VRTRasterBand dataType="CFloat32" band="{index}">
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
    </VRTRasterBand>\n'''

class Stack(object):

    def __init__(self, slcDir=None):
        self.slcDir = slcDir
        self.bbox = None
    def configure(self, outDir):
        if not os.path.exists(outDir):
            os.makedirs(outDir)

        self.outSlcVrtDir = os.path.join(outDir, "slcs")
        self.outStackVrtDir = os.path.join(outDir, "stack")

        if not os.path.exists(self.outSlcVrtDir):
            os.makedirs(self.outSlcVrtDir)
        
        if not os.path.exists(self.outStackVrtDir):
            os.makedirs(self.outStackVrtDir)
        
    def gatherSLCs(self):
        
        self.slcList = glob.glob(os.path.join(self.slcDir, '*/*.slc'))
        if len(self.slcList) == 0:
           self.slcList = glob.glob(os.path.join(self.slcDir, '*/*.slc.full'))

        print('Number of SLCs discovered: ', len(self.slcList))
        print('We assume that the SLCs and the ann files are sorted in the same order')

        self.slcList.sort()
        self.size = len(self.slcList)
        self.applyBbox = [True]*self.size

    def getSize(self):
        self.size = len(self.slcList)

    def getDates(self):

        self.dateList = []
        for ind, slc in enumerate(self.slcList):
            dd = os.path.basename(os.path.dirname(slc))
            self.dateList.append(dd)

    def writeStackVRT(self):
        # modified this function from stripmap2vrt.py
        dates = []
        data = [] 
        for ind, slc in enumerate(self.slcList):

            metadata = {}
            width = None
            height = None
            path = None

            print("*****")
            print(slc)

            ## create xml file if missing
            #if not os.path.isfile(slc+'.xml'):
            #    cmd = 'gdal2isce_xml.py -i '+slc
            #    print(cmd)
            #    os.system(cmd)

            ## fix potential filepath in xml file if dir has been moved.
            #img = IML.loadImage(slc)[0]
            #img.filename = slc
            #img.setAccessMode('READ')
            #img.renderHdr()

            ds = gdal.Open(slc + '.vrt', gdal.GA_ReadOnly)
            width = ds.RasterXSize
            height = ds.RasterYSize
            ds = None


            metadata['WAVELENGTH'] = 0.03
            metadata['ACQUISITION_TIME'] = os.path.basename(os.path.dirname(slc))

            path = os.path.abspath(slc)

            tag = metadata['ACQUISITION_TIME']
            
#        outname =  datetime.datetime.strptime(tag.upper(), '%d-%b-%Y %H:%M:%S UTC').strftime('%Y%m%d')

            outname = metadata['ACQUISITION_TIME']
            with open( os.path.join(self.outSlcVrtDir, '{0}.vrt'.format(outname)) , 'w') as fid:
                fid.write( vrttmpl.format(width=width,
                                     height=height,
                                     PATH=path,
                                     linewidth=8*width))


            data.append(metadata)
            dates.append(outname)

        '''
        # setting up a subset of the stack
        ymin, ymax, xmin, xmax = [0 , height, 0 , width]
        if self.bbox:
            ymin, ymax, xmin, xmax = self.bbox

        xsize = xmax - xmin
        ysize = ymax - ymin
        '''

        self.stackVRT = os.path.join(self.outStackVrtDir, 'stack.vrt')
        print("writing ", self.stackVRT)
        with open( self.stackVRT, 'w') as fid:
        #fid.write( '<VRTDataset rasterXSize="{width}" rasterYSize="{height}">\n'.format(width=width, height=height))
            width, height, xmin, ymin, xsize, ysize = self.get_x_y_offsets(-1)
            fid.write( '<VRTDataset rasterXSize="{xsize}" rasterYSize="{ysize}">\n'.format(xsize=xsize, ysize=ysize))

            for ind, (date, meta) in enumerate( zip(dates, data)):
                width, height, xmin, ymin, xsize, ysize = self.get_x_y_offsets(ind)
                outstr = vrttmpl2.format(width=width, height=height,
                                xmin=xmin, ymin=ymin,
                                xsize=xsize, ysize=ysize,
                                date=date, acq=meta['ACQUISITION_TIME'],
                                wvl = meta['WAVELENGTH'], index=ind+1,
                                path = os.path.abspath( os.path.join(self.outSlcVrtDir, date+'.vrt')))
                fid.write(outstr)

            fid.write('</VRTDataset>')
        

    def get_x_y_offsets(self, ind):
        # setting up a subset of the stack
        slc = self.slcList[ind]
        ds = gdal.Open(slc + '.vrt', gdal.GA_ReadOnly)
        width = ds.RasterXSize
        height = ds.RasterYSize
        ds = None

        ymin, ymax, xmin, xmax = [0 , height, 0 , width]
        print(ind, self.applyBbox[ind])
        if (self.bbox): 
            if self.applyBbox[ind]:
                ymin, ymax, xmin, xmax = self.bbox

        xsize = xmax - xmin
        ysize = ymax - ymin
         
        return width, height, xmin, ymin, xsize, ysize

class MiniStack(Stack):

    #def __init__(self, slcList=None):
    #    self.slcList = slcList

    def updateMiniStack(self, compressedSlcDir):
        compSlcList = glob.glob(os.path.join(compressedSlcDir, '*/*.slc'))
        applyBbox = [False]*len(compSlcList)
        self.slcList =  compSlcList + self.slcList

        self.applyBbox = applyBbox + self.applyBbox  # Note: order matters here. 
    

