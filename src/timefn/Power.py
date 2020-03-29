#!/usr/bin/env python3

############################################################
# Copyright 2013, by the California Institute of Technology#
# Author : Piyush Agram
############################################################

from .CenteredBasisFn import CenteredBasisFn
import numpy as np 

class Power(CenteredBasisFn):
    '''Class to represent power function.

    .. math ::

        (\frac{t_{norm}}{\tau})^{exp}
    '''

    fnname = 'Power'

    def __init__(self, tau=1, exp=None, **kwargs):
        '''Constructor for the power function

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
            Time scale for the power function
        exp: str (or) flat
            Exponent for the power function
        side: str
            Alternate interface to tmin, tmax
        '''
       
        ####Initialize the base class
        CenteredBasisFn.__init__(self, **kwargs)
        
        try:
            self.exponent = float(exp)
        except:
            raise TypeError('Exponent {0} cannot be interpreted as float'.format(exp))

        try:
            self.tau = float(tau)
        except:
            raise TypeError('Tau {0} cannot be interpreted as float'.format(tau))

        if self.tau == 0:
            raise ValueError('Normalizing factor Tau cannot be zero.')

        self.fixedExp = False

    
    def computeRaw(self, tinput):
        '''Evaluation of the power function.
        '''
        return  np.power(self.normalizedTime(tinput)/self.tau, self.exponent)

    def kwrepr(self):
        '''Keyword representation.
        '''
        outstr = CenteredBasisFn.kwrepr(self)
        if outstr:
            outstr += ','

        if not self.fixedExp:
            outstr += 'exp=' + str(self.exponent)+','

        if self.tau != 1:
            outstr += 'tau=' + str(self.tau)

        return outstr


class Linear(Power):
    '''Class for linear function.

    This is same as the Power function, with exp set to 1.0.
    '''

    fnname = 'Linear'

    def __init__(self, **kwargs):
        '''Constructor for Linear function.

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
            Time scale for the power function
        side: str
            Alternate interface to tmin, tmax
        '''

        Power.__init__(self, exp=1.0, **kwargs)
        self.fixedExp = True


class Quadratic(Power):
    '''Class for quadratic function.

    This is the same as the Power function, with exp set to 2.0.
    '''

    fnname = 'Quadratic'

    def __init__(self, **kwargs):
        '''Constructor for Quadratic function.

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
            Time scale for the power function
        side: str
            Alternate interface to tmin, tmax
        '''
        Power.__init__(self, exp=2.0, **kwargs)
        self.fixedExp = True


class Cubic(Power):
    '''Class for cubic function.

    This is same as the Power function, with exp set to 3.0.
    '''

    fnname = 'Cubic'

    def __init__(self, **kwargs):
        '''Constructor for Cubic function.

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
            Time scale for the power function
        side: str
            Alternate interface to tmin, tmax
        '''
        Power.__init__(self, exp=3.0, **kwargs)
        self.fixedExp = True


class Quartic(Power):
    '''
    Class for quartic function.

    This is the same as Power function with exp set to 4.0.
    '''

    fnname = 'Quartic'

    def __init__(self, **kwargs):
        '''Constructor for Quartic function.

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
            Time scale for the power function
        side: str
            Alternate interface to tmin, tmax
        '''
        Power.__init__(self, exp=4.0, **kwargs)
        self.fixedExp = True
