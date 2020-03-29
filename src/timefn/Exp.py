#!/usr/bin/env python3

############################################################
# Copyright 2013, by the California Institute of Technology#
# Author : Piyush Agram
############################################################


from .CenteredBasisFn import CenteredBasisFn
import numpy as np 

class Exp(CenteredBasisFn):
    '''Class to represent exp function (base e).

    .. math ::

        e^(\frac{t_{norm}}{\tau}

    '''

    fnname = 'Exp'

    def __init__(self, tau=None, **kwargs):
        '''Constructor for the exp function.

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
            Time scale for the exp function (keyword)
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
        '''Computes the exponential function.
        '''
        return  np.exp( - self.efactor * self.normalizedTime(tinput) / self.tau)


    def kwrepr(self):
        '''Keyword representation as a string.
        '''
        outstr = CenteredBasisFn.kwrepr(self)
        if outstr:
            outstr += ','

        return  outstr+ 'tau=' + str(self.tau)


class Exp10(Exp):
    '''Class to represent exp function (base 10).

    .. math ::

        10^(\frac{t_{norm}}{\tau})

    '''

    fnname = 'Exp10'

    def __init__(self, **kwargs):
        '''Constructor for the exp10 function.

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
            Time scale for the exp10 function (keyword)
        side: str
            Alternate interface to tmin, tmax (keyword)

        '''
       
        ####Initialize base class
        Exp.__init__(self, **kwargs)
        
        self.efactor = np.log(10.0)


class Exp2(Exp):
    '''Class to represent exp function (base 2).

    .. math :: 

        2^(\frac{t_{norm}}{\tau})

    '''
    
    fnname = 'Exp2'

    def __init__(self, **kwargs):
        '''Constructor for the exp2 function.

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
            Time scale for the exp2 function (keyword)
        side: str
            Alternate interface to tmin, tmax (keyword)
        '''

        ####Initialize base class
        Exp.__init__(self, **kwargs)

        self.efactor = np.log(2.0)
