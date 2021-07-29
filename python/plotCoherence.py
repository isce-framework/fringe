#!/usr/bin/env python3

from osgeo import gdal
import argparse
import numpy as np
import matplotlib.pyplot as plt

def cmdLineParse():
    '''
    Command line parse.
    '''

    parser = argparse.ArgumentParser(description='Test MLE for single pixel')
    parser.add_argument('-i', dest='inputDS', type=str, required=True,
            help='Input Stack VRT')
    parser.add_argument('-w', dest='wtsDS', type=str, required=True,
            help='Weights dataset')
    parser.add_argument('-y', dest='line', type=int, required=True,
            help='Line number of pixel of interest')
    parser.add_argument('-x', dest='pixel', type=int, required=True,
            help='Pixel number of pixel of interest')
    return parser.parse_args()

class Dummy(object):
    pass

class BitMask:
    def __init__(self, Ny, Nx):
        self.Ny = Ny
        self.Nx = Nx

    def getbit(self, mask, ii,jj):
        flat = (ii+self.Ny) * (2*self.Nx+1) + jj + self.Nx
        num = flat//8
        bit = flat % 8

        return ((mask[num] >> bit) & 1) 


def unpack(seq):
    bits = []
    for x in seq:
        for kk in range(8):
            bits.append((x>>kk) & 1)
    return bits


def loadData(inps):
    '''
    Load relevant data for a pixel.
    '''

    ds = gdal.Open(inps.wtsDS, gdal.GA_ReadOnly)
    Nx = int(ds.GetMetadataItem('HALFWINDOWX', 'ENVI'))
    Ny = int(ds.GetMetadataItem('HALFWINDOWY', 'ENVI'))
    width = ds.RasterXSize
    lgth = ds.RasterYSize
    bands = ds.RasterCount
    ds = None

    fid = open(inps.wtsDS, 'rb')
    fid.seek( (inps.line * width + inps.pixel) * bands * 4)
    mask = fid.read(bands*4)
    fid.close()

    #bits = unpack(mask)
    npix = (2*Ny+1) * (2*Nx+1)
    #bits = bits[:npix]
    #bitmask=(np.array(bits) == 1)
    #count = np.sum(bits)
    masker = BitMask(Ny,Nx)

    bitmask = np.zeros(npix, dtype=bool)
    count = 0
    ind = 0
    for ii in range(-Ny,Ny+1):
        for jj in range(-Nx,Nx+1):
            flag = masker.getbit(mask, ii, jj)
            print("ii = {0}, jj={1}, flag = {2}".format(ii,jj,flag))
            bitmask[ind] = (flag == 1)
            if flag:
                count += 1
            ind += 1


    print("Line: ", inps.line)
    print("Pixel: ", inps.pixel)
    print("Count: ", count)

    xoff = inps.pixel - Nx
    yoff = inps.line - Ny
  
    ds = gdal.Open(inps.inputDS, gdal.GA_ReadOnly)
    nslc = ds.RasterCount

    print("Stack size: ", nslc)
    data = np.zeros((count, nslc), dtype=np.complex64)
    for slc in range(nslc):
        indata = ds.GetRasterBand(slc+1).ReadAsArray(xoff, yoff, 2*Nx+1, 2*Ny+1).flatten()
        data[:,slc] = indata[bitmask]

    ds = None
    return data

def covariance(C1,C2):
    A1 = np.sum(np.abs(C1)**2)
    A2 = np.sum(np.abs(C2)**2)
    cov = np.sum(C1*np.conjugate(C2))/(np.sqrt(A1)*np.sqrt(A2))
    return cov


def computeCovar(stack):
    numberOfSamples , numberOfSlc = stack.shape
    cov_mat = np.zeros((numberOfSlc , numberOfSlc), dtype=np.complex64)
    for ti in range(numberOfSlc):
        for  tj in range(ti+1, numberOfSlc):
            cov = covariance(stack[:,ti], stack[:,tj])
            cov_mat[ti,tj] = cov
            cov_mat[tj,ti] = np.conjugate(cov)
        cov_mat[ti,ti] = 1.0

    return cov_mat

if __name__ == '__main__':
    '''
    Main driver.
    '''
    
    ###Parse command line
    inps = cmdLineParse()

    ###Read in relevant data
    data = loadData(inps)

    ###Compute covariance
    cov = computeCovar(data)

    coh = np.abs(cov)
    
    fig = plt.figure()
    ax = fig.add_subplot(1,1,1)
    im = ax.imshow(np.abs(cov), vmin = 0.0, vmax=1 )
    fig.colorbar(im, ax=ax)

    plt.show()

    


