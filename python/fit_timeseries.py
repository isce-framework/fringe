#!/usr/bin/env python3

# Author: Heresh Fattahi

import os
import glob
import argparse
import numpy as np
from osgeo import gdal
import time
import datetime
from lmfit import Minimizer, Parameters, report_fit
import multiprocessing as mp


PI = np.pi
wvl = 0.05546576

params = Parameters()
params.add('lnTau', value=0, min = -10, max = 1)
params.add('P', value=0)
params.add('C', value=0)
params.add('offset', value = 0)
params.add('rate',value = 0)


def cmdLineParser():
    '''
    Command line parser.
    '''

    parser = argparse.ArgumentParser(description = 'non linear fit to a time-series.',
                formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument('-t', '--time_series_file', type=str, dest='timeSeries',
            required=True, help='input file for time-series')

    parser.add_argument('-n', '--nproc', type=int, dest='nproc', default = 8,
            help='number of processors')

    return parser.parse_args()

def get_time(timeSeries_file):

    ds = gdal.Open(timeSeries_file, gdal.GA_ReadOnly)

    nslc = ds.RasterCount
    print("Stack size: ", nslc)
    t = np.zeros(nslc)
    for slc in range(nslc):
        m = ds.GetRasterBand(slc+1).GetMetadata("slc")
        t[slc] = m['AcquisitionTime']
    ds = None
    return t

def get_data_chunk(ds , nslc, width, y, ysize):
    dataCube = ds.ReadAsArray(0, y, width, ysize)
    return  dataCube*wvl/4.0/PI

def get_data(ds , data, nslc, x, y):

    for slc in range(nslc):
        indata = ds.GetRasterBand(slc+1).ReadAsArray(x, y, 1, 1)
        data[slc] = indata
    return  data*wvl/4.0/PI


def residual(params, t, data, t_eq):
    offset = params['offset']
    rate = params['rate']
    lnTau = params['lnTau']
    C = params['C']
    P = params['P']
    post = np.log10(1+(t-t_eq)/np.exp(lnTau))
    post[np.isnan(post)]=0
    post[post==-np.inf]=0
    post[post==np.inf]=0
    H = 1 * ((t-t_eq) > 0)
    model = offset + rate*t + H*(C + P*post)
    #model = offset + rate*t + heaviside(t,t_eq)*(C + P*post)
    return model - data

#@jit('f8,f8')
def heaviside(t,t_eq):
    x = t - t_eq
    return 1 * (x > 0)


class fit_data_multiproc(mp.Process):
    
    def __init__(self, par):
        self.par = par
        mp.Process.__init__(self)

    def run(self):

        line_inds = self.par.line_inds
        nlines = len(line_inds)
        pixels = self.par.pixels

        for q in range(nlines):
            for jj in range(pixels):
                minner = Minimizer(residual, self.par.params, fcn_args=(self.par.t, self.par.data[:,q,jj], self.par.t_EQ))
                result = minner.minimize()

                ii = line_inds[q]
                self.par.Co[ii, jj] = result.params['C'].value
                self.par.P[ii, jj] = result.params['P'].value
                self.par.Rate[ii, jj] = result.params['rate'].value
                self.par.Off[ii, jj] = result.params['offset'].value
                self.par.lnTau[ii, jj] = result.params['lnTau'].value 
            

class dummy:
    pass
      


def fit_data(inps, params, t_EQ):
 
    nproc = inps.nproc 
    #t = get_time("/u/k-data/fattahi/Kurdistan/TimeSeries/72/Sequential/TimeSeries_Post/timeSeries.vrt")
    t = get_time(inps.timeSeries)
    ds = gdal.Open(inps.timeSeries, gdal.GA_ReadOnly)
    nslc = ds.RasterCount
    lines = ds.RasterYSize
    pixels = ds.RasterXSize

    
    dsC = np.memmap("Co_seismic.bin", dtype=np.float32, mode='w+', shape=(lines,pixels))
    dsP = np.memmap("Post_seismic.bin", dtype=np.float32, mode='w+', shape=(lines,pixels))
    dsRate = np.memmap("rate.bin", dtype=np.float32, mode='w+', shape=(lines,pixels))
    dsOff = np.memmap("offset.bin", dtype=np.float32, mode='w+', shape=(lines,pixels))
    dsTau = np.memmap("lnTau.bin", dtype=np.float32, mode='w+', shape=(lines,pixels))
    

    ind = [ii for ii in range(0, lines, 128)]
    if ind[-1] != lines:
       ind.append(lines)
    print(ind)


    ind_diff = np.diff(ind)
    ind = ind[0:-1]
    print(ind)
    print(ind_diff)

    pinds = np.int_(np.linspace(0,pixels,num=nproc+1))
    par = dummy()
    par.pixels = pixels
    par.params = params
    par.t_EQ = t_EQ
    par.t = t 
    for ii,line in enumerate(ind):
        print("at line : ", line)
        threads = []
        line_inds = np.int_(np.linspace(0, ind_diff[ii], num=nproc+1))
        dataChunk = get_data_chunk(ds , nslc, pixels, line, int(ind_diff[ii]))

        # shared memory objects for estimated parameters (1 chunk in azimuth x width in range)
        tempCo = mp.Array('d',int(ind_diff[ii])*pixels)
        tempP = mp.Array('d',int(ind_diff[ii])*pixels)
        tempRate = mp.Array('d',int(ind_diff[ii])*pixels)
        tempOff = mp.Array('d',int(ind_diff[ii])*pixels)
        templnTau = mp.Array('d',int(ind_diff[ii])*pixels)

        par.Co = np.reshape(np.frombuffer(tempCo.get_obj()),(ind_diff[ii], pixels))
        par.P  = np.reshape(np.frombuffer(tempP.get_obj()),(ind_diff[ii], pixels))
        par.Rate = np.reshape(np.frombuffer(tempRate.get_obj()),(ind_diff[ii], pixels))
        par.Off = np.reshape(np.frombuffer(tempOff.get_obj()),(ind_diff[ii], pixels))
        par.lnTau = np.reshape(np.frombuffer(templnTau.get_obj()),(ind_diff[ii], pixels))

        ts = time.time()
        for q in range(nproc):
                inds = np.arange(line_inds[q], line_inds[q+1])
                par.line_inds = inds
                par.data = dataChunk[:,inds,:]
                threads.append(fit_data_multiproc(par))
                threads[q].start()

        for thrd in threads:
                thrd.join()
        print(time.time() - ts , " sec")
         
        dsC[line:line+ind_diff[ii] , :] = par.Co
        dsP[line:line+ind_diff[ii] , :] = par.P
        dsRate[line:line+ind_diff[ii] , :] = par.Rate
        dsOff[line:line+ind_diff[ii] , :] = par.Off
        dsTau[line:line+ind_diff[ii] , :] = par.lnTau
        
    dsC = None
    dsP = None
    dsRate = None
    dsOff = None
    dsTau = None
    ds = None          

if __name__ == '__main__':
    '''
    Main driver.
    '''
    #*************************************************************#
    # read the input options and unwrap
    inps = cmdLineParser()
    
    params = Parameters()
    params.add('lnTau', value=0, min = -10, max = 1)
    params.add('P', value=0)
    params.add('C', value=0)
    params.add('offset', value = 0)
    params.add('rate',value = 0)

    teq="2017-11-12 18:18:17 UTC"
    deq = datetime.datetime.strptime(teq.upper(), '%Y-%m-%d %H:%M:%S UTC')
    day_of_year = deq.timetuple().tm_yday
    t_EQ = float(deq.year)+float(day_of_year-1)/365.25


    fit_data(inps, params, t_EQ)




