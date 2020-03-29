#!/usr/bin/env python3

############################################################
# Copyright 2013, by the California Institute of Technology#
# Author : Piyush Agram
############################################################

from .CenteredBasisFn import CenteredBasisFn
import numpy as np 
from scipy.misc import factorial, comb

class BSpline(CenteredBasisFn):
    '''Class to represent a Bspline.

    Attributes
    ----------
    scale : float
        Time scale of the BSpline
    order : int
        Order of Bspline (default = 2)

    Notes
    -----

    For description : http://mathworld.wolfram.com/B-Spline.html

    tref is of more importance than tmin, tmax for BSplines.
    For now, avoid using the tmin, tmax keywords.
    '''

    fnname = 'BSpline'

    def __init__(self, scale=None, order=2, **kwargs):
        '''Constructor for Bspline function.

        Parameters
        ----------
        tref: str (or) datetime.datetime
            Reference epoch (time zero) (keyword)
        units: str
            Units to interpret tau (keyword)
        scale: str (or) float
            Time scale for the spline function
        order: str (or) int
            Order of the spline function

        '''        
        CenteredBasisFn.__init__(self, **kwargs)
        
        try:
            self.scale = float(scale)
        except:
            raise TypeError('Scale {0} cannot be interpreted as float'.format(scale))

        if self.scale == 0 :
            raise ValueError('Time scale of a BSpline cannot be zero.')


        try:
            self.order = int(order)
        except:
            raise TypeError('Order {0} cannot be interpreted as int'.format(order))

        if self.order < 2:
            raise ValueError('Order of a BSpline cannot be less than 2.')

    
    def computeRaw(self, tinput):
        '''Bspline computation.
        '''
        tnorm = self.normalizedTime(tinput) / self.scale
        x = tnorm + self.order + 1
        b = 0.0
        for k in range(self.order+2):
            m  = x-k-(self.order+1)/2
            up = np.power(m,self.order)
            b += ((-1)**k)*comb(self.order+1,k)*up*(m>=0)

        b = b / (1.0*factorial(self.order))
        return b

    def kwrepr(self):
        '''Keyword representation.
        '''
        outstr = CenteredBasisFn.kwrepr(self)
        if outstr:
            outstr += ','

        outstr += 'scale=' + str(self.scale)

        if self.order != 2:
            outstr += ',order=' + str(self.order) 

        return outstr

class ISpline(CenteredBasisFn):
    '''Class to represent a ISpline.

    Attributes
    ----------
    scale : float
        Time scale of the ISpline
    order : int
        Order of ISpline (default = 2)

    Notes
    -----

    ISpline is an integral of a BSpline.

    tref is of more importance than tmin, tmax for ISplines.
    For now, avoid using the tmin, tmax keywords.
    '''

    fnname = 'ISpline'

    def __init__(self, scale=None, order=2, **kwargs):
        '''Constructor for ISpline function.

        Parameters
        ----------
        tref: str (or) datetime.datetime
            Reference epoch (time zero) (keyword)
        units: str
            Units to interpret tau (keyword)
        scale: str (or) float
            Time scale for the spline function
        order: int
            Order of the spline function

        '''          
        CenteredBasisFn.__init__(self, **kwargs)
        
        try:
            self.scale = float(scale)
        except:
            raise TypeError('Scale {0} cannot be interpreted as float'.format(scale))

        if self.scale == 0 :
            raise ValueError('Time scale of a ISpline cannot be zero.')

        try:
            self.order = int(order)
        except:
            raise TypeError('Order {0} cannot be interpreted as int'.format(order))

        if self.order < 2:
            raise ValueError('Order of a ISpline cannot be less than 2.')

    
    def computeRaw(self, tinput):
        '''Ispline computation.
        '''
        tnorm = self.normalizedTime(tinput) / self.scale
        x = tnorm + self.order + 1
        b = 0.0
        for k in range(self.order+2):
            m = x-k-(self.order+1)/2
            up = np.power(m,self.order+1)
            b += ((-1)**k)*comb(self.order+1,k)*up*(m>=0)

        b = b * self.scale /( (self.order+1.0) * factorial(self.order))
        return b

    def kwrepr(self):
        '''Keyword representation.
        '''
        outstr = CenteredBasisFn.kwrepr(self)

        if outstr:
            outstr += ','

        outstr += 'scale=' + str(self.scale)
        if self.order != 2:
            outstr += ',order=' + str(self.order) 

        return outstr


class PBSpline(BSpline):
    '''Class to represent a Periodic Bspline.

    Attributes
    ----------
    scale : float
        Time scale of the Periodic BSpline
    order : int
        Order of Bspline (default = 2)
    period : float
        Time period for the Periodic BSpline
    offset : float
        Offset within a time period.

    Notes
    -----
    tref is of more importance than tmin, tmax for BSplines.
    For now, avoid using the tmin, tmax keywords.
    '''

    fnname = 'PBSpline'

    def __init__(self, period=None, offset=None, **kwargs):
        '''Constructor for Periodic Bspline function.

        Parameters
        ----------
        tref: str (or) datetime.datetime
            Reference epoch (time zero) (keyword)
        units: str
            Units to interpret tau (keyword)
        scale: str (or) float
            Time scale for the spline function
        period: str (or) float
            Time period for the spline function
        offset: str (or) float
            Offset within a time period
        order: str (or) int
            Order of the spline function

        '''        
        BSpline.__init__(self, **kwargs)
        
        try:
            self.period = float(period)
        except:
            raise TypeError('Period {0} cannot be interpreted as float.'.format(period))

        if self.period == 0:
            raise ValueError('Period of a PBSpline cannot be zero.')

        try:
            self.offset = float(offset)
        except:
            raise TypeError('Offset {0} cannot be interpreted as float.'.format(offset))

        if self.offset >= self.period:
            raise ValueError('Offset {0} cannot exceed the time period {1} for Periodic BSplines'.format(self.offset, self.period))


    def normalizedTime(self, tinput):
        '''Extending the CenteredBasisFn method to account for offset and period.
        '''

        ###Usual normalized time
        tnorm = BSpline.normalizedTime(self, tinput)

        ###Adjust for offset and time period
        tnorm = tnorm - self.offset
        tnorm = tnorm - self.period * np.round(tnorm/self.period)

        return tnorm
    
    def kwrepr(self):
        '''Keyword representation.
        '''

        outstr = BSpline.kwrepr(self)
        outstr += ',period=' + str(self.period)

        return outstr

