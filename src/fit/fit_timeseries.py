#!/usr/bin/env python3
# ----------------------------------------------------------------------------
# Copyright (c) 2018-, California Institute of Technology ("Caltech"). 
# U.S. Government sponsorship acknowledged.
# All rights reserved.
# 
# Author(s): Heresh Fattahi
# ----------------------------------------------------------------------------

import os
import argparse
import datetime
import time
import gdal
import numpy as np
import fitlib as Fit
from TimeFunction import TimeFunction

def cmdLineParser():
    '''
    Command line parser.
    '''

    parser = argparse.ArgumentParser(description = 'Fitting a time-series with a mathematical model',
                formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('-i', '--input', type=str, dest='inputDS',
            required=True, help='Input dataset (A single VRT, whose bands are different acquisition dates)')
    parser.add_argument('-o', '--output', type=str, dest='outputDS',
            required=True, help='Output dataset')
    parser.add_argument('-t', '--earthquakesTimes', type=str, nargs = '+', dest='earthquakesTimes',
            default=None, help='A list of the time of earthquake YYYY-MM-DD')
    parser.add_argument('-T', '--tau', type=float, nargs='+', dest='tau',
            default=None, help='A list of Tau values in days (relaxation time of the post-seismic signal)')
    parser.add_argument('-p', '--periods', type=float, nargs='+', dest='periods',
            default=None, help='A list of the periods of periodic terms in days (e.g., for annual and semi-annual terms: -p 365 182.5)') 
    parser.add_argument('-P', '--polynomial_order', type=int, dest='polynomial_order',
            default=1, help='order of the polynomial (default: 1)')
    parser.add_argument('-b', '--blockSize', type=int, dest='blockSize',
            default=25, help='number of lines for each processing block of data')

    return parser.parse_args()   
 
def runFit(inps, G):

    fitObj = Fit.PyFit()
    fitObj.set_blocksize(inps.blockSize)
    fitObj.set_outputDataset(inps.outputDS)
    fitObj.open_inputDataset(inps.inputDS)
    fitObj.set_G(G)
    fitObj.set_Npar(G.shape[1])
    fitObj.fit_timeseries()

def extract_time(inputDS):
    ds = gdal.Open(inputDS, gdal.GA_ReadOnly)
    nBands = ds.RasterCount
    time_array=[]
    for bnd in range(nBands):
        dd = ds.GetRasterBand(bnd+1).GetMetadata("slc")['Date']
        dd = datetime.datetime(*time.strptime(dd, "%Y%m%d")[0:6])
        time_array.append(dd) 

    return np.array(time_array)

def getG(time_array, t_events, tau_events, tref=None, poly_order=None, periods=[]):

    if tref is None:
       tref = time_array[0]

    # instatiate a TimeFunction object
    tfn = TimeFunction(time_array=time_array, units='days') 
    tfn.configure()

    if poly_order:
        tfn.add_polynomial(tref, poly_order) 

    # add step and logarithmic functions for multiple events
    for ii, t in enumerate(t_events):
        tfn.add_heaviside(t)
        if tau_events[ii]>0:
            tfn.add_post_seismic(t, tau_events[ii])
        else:
            print("Tau can not be zero. Ignoring post-seismic for {0}".format(t))

    # add periodic terms
    if periods:
        for p in periods:
            tfn.add_periodic(time_array[0], p)

    G, timeFnStr = tfn.getG()
 
    return G, timeFnStr

if __name__ == '__main__':
    '''
    Main driver.
    '''
    inps = cmdLineParser()
    inps.inputDS = os.path.abspath(inps.inputDS)
    inps.outputDS = os.path.abspath(inps.outputDS)

    # extract the acquisition times from the input vrt file
    time_array = extract_time(inps.inputDS)

    # create a list of event times
    earthquake_times = []
    if inps.earthquakesTimes is not None:
        for eq in inps.earthquakesTimes:
            eqt = datetime.datetime(*time.strptime(eq, "%Y-%m-%d")[0:6])
            earthquake_times.append(eqt)

    # Build the G matrix for the model ( d = Gm)
    G, timeFnStr = getG(time_array, earthquake_times, inps.tau, 
                         periods=inps.periods, 
                         poly_order = inps.polynomial_order)

    print("****************************")
    print("Requested time function has the following form: ")
    print("")
    print(timeFnStr)
    print("")
    print("****************************")

    # estimating the model parameters
    runFit(inps, G)

