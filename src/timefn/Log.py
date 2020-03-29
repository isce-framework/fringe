#!/usr/bin/env python3

############################################################
# Copyright 2013, by the California Institute of Technology#
# Author : Piyush Agram
############################################################

from .CenteredBasisFn import CenteredBasisFn
import numpy as np 

class Log(CenteredBasisFn):
    '''Class to represent log function (base e).

    .. math ::

        ln(|\frac{t_{norm}}{\tau}|)
    '''

    fnname = 'Log'

    def __init__(self, tau=None, **kwargs):
        '''Constructor for the natural logarithm function.

        Parameters
        ----------
        tmin: str (or) datetime.datetime
            Left edge of active period (keyword)
        tmax: str (or) datetime.datetime
            Right edge of active period (keyword)
        tref: str (or) datetime.datetime
            Reference epoch (time zero) (keyword)
        units: str
            Units to interpret tau (keyword) 
        tau: str (or) float
            Time scale for the log function (keyword)
        side: str
            Alternate interface to tmin, tmax (keyword)
        '''
        
        CenteredBasisFn.__init__(self, **kwargs)
        
        try:
            self.tau = float(tau)
        except:
            raise TypeError('Tau {0} cannot be interpreted as float'.format(tau))

        
        if self.tau == 0:
            raise ValueError('Normalizing factor Tau cannot be zero.')

        self.efactor = 1.0
    
    def computeRaw(self, tinput):
        '''Computes the log function.

        Notes
        -----

        We use absolute value here to avoid issues with computing log of negative values.
        '''

        return  np.log( np.abs(self.normalizedTime(tinput) / self.tau)) / self.efactor

    def kwrepr(self):
        '''Keyword representation.
        '''
        outstr = CenteredBasisFn.kwrepr(self)
        if outstr:
            outstr += ','
        return  outstr + 'tau=' + str(self.tau)

class Log10(Log):
    '''Class to represent log function (base 10).

    .. math ::

        log_{10}(|\frac{t_{norm}}{\tau}|)
    '''

    fnname = 'Log10'

    def __init__(self, **kwargs):
        '''Constructor for the Log10 function.

        Parameters
        ----------
        tmin: str (or) datetime.datetime
            Left edge of active period (keyword)
        tmax: str (or) datetime.datetime
            Right edge of active period (keyword)
        tref: str (or) datetime.datetime
            Reference epoch (time zero) (keyword)
        units: str
            Units to interpret tau (keyword)
        tau: str (or) float
            Time scale for the log function (keyword)
        side: str
            Alternate interface to tmin, tmax (keyword)
        '''
        
        Log.__init__(self, **kwargs)
        
        self.efactor = np.log(10.0)


class Log2(Log):
    '''Class to represent log function (base 2).

    .. math ::

        log_{2}(|\frac{t_{norm}}{\tau}|)
    '''

    fnname = 'Log2'

    def __init__(self, **kwargs):
        '''Constructor for Log2 function

        Parameters
        ----------
        tmin: str (or) datetime.datetime
            Left edge of active period (keyword)
        tmax: str (or) datetime.datetime
            Right edge of active period (keyword)
        tref: str (or) datetime.datetime
            Reference epoch (time zero) (keyword)
        units: str
            Units to interpret tau (keyword)
        tau: str (or) float
            Time scale for the log function (keyword)
        side: str
            Alternate interface to tmin, tmax (keyword)
        '''

        Log.__init__(self, **kwargs)

        self.efactor = np.log(2.0)


class GeoLog(Log):
    '''Class to use for log functions (geophysical models).

    .. math ::

        ln( 1.0 + | \frac{t_{norm}}{\tau}|)
    '''

    fnname = 'GeoLog'
    
    def __init__(self, **kwargs):
        '''Constructor for GeoLog function.

        Parameters
        ----------
        tmin: str (or) datetime.datetime
            Left edge of active period (keyword)
        tmax: str (or) datetime.datetime
            Right edge of active period (keyword)
        tref: str (or) datetime.datetime
            Reference epoch (time zero) (keyword)
        units: str
            Units to interpret tau (keyword)
        tau: str (or) float
            Time scale for the geolog function (keyword)
        side: str
            Alternate interface to tmin, tmax (keyword)
        '''

        kwargs['side'] = 'right'  ####Explicitly set to right
        Log.__init__(self, **kwargs)

    def computeRaw(self, tinput):
        '''Compute the raw function.
        '''
        return np.log(1.0 + np.abs(self.normalizedTime(tinput)/self.tau))
