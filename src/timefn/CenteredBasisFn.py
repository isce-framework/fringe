#!/usr/bin/env python3

############################################################
# Copyright 2013, by the California Institute of Technology#
# Author : Piyush Agram
############################################################

from .BasisFn import BasisFn
from . import datetimeUtils as DU
import numpy as np

try: 
    from pint import UnitRegistry
except ImportError:
    raise ImportError('Python pint library not available. Please install pint before proceeding.')

ureg = UnitRegistry()

class CenteredBasisFn(BasisFn):
    ''' Base class for timefn that needs a reference epoch in time.

    Attributes
    ----------
    tref : datetime.datetime
        Reference epoch for evluation of the function.
    units : str
        Time scale to be used for interpreting model parameters.
        

    Notes
    -----

    CenteredBasisFn is derived from BasisFn.
    

    '''

    fnname = 'CenteredBasis'

    def __init__(self, tref=None, units='secs', side='both', **kwargs):
        ''' Constructor of the CenteredBasisFn class.

        Parameters
        ----------

        tmin: str (or) datetime.datetime
            Left edge of active period of function (keyword)
        tmax: str (or) datetime.datetime
            Right edge of active period of function (keyword)
        tref: str (or) datetime.datetime
            Reference epoch for evaluation of function (keyword)
        units: str
            String representing the scale. Python pint is used to translate units to SI units internally for all computation.
        side: str
            Can be left, right or both - referring to active period w.r.t reference epoch.
        '''

        ####Initialize the base class
        BasisFn.__init__(self, **kwargs)

        self.tref = DU.interpretAsDatetime(tref)

        #####Convert to SI units 
        try:
            self.norm = float(ureg[units].to_base_units().magnitude)
        except:
            raise ValueError('Undefined units {0} provided. Cannot convert to SI units'.format(units))

        self.units = units
        
        if side not in ['left', 'right', 'both']:
            raise Exception('Side {0} not valid. Should be left / right / both. '.format(side))
        elif side == 'left':
            self.tmax = tref
        elif side == 'right':
            self.tmin = tref

        ####Use numpy vectorize to handle - scalar or an array
        self.normalizedTime = np.vectorize(self.normalizedTimeOne)

    def normalizedTimeOne(self, tinput):
        '''Return normalized time. Difference from tref and scaled by units.

        Parameters
        ----------
        tinput : datetime.datetime
            Returns normalized time for computing the functions.

        .. math:: 
        
            t_{norm} = (t - t_{ref})/units

        '''

        return (tinput-self.tref).total_seconds() / self.norm
#        return  np.array([(x-self.tref).total_seconds() for x in tinput]) / self.norm


    def kwrepr(self):
        '''String representation of function keywords.

        Note
        ----
        In addition to tmin, tmax, we now include tref and units.

        '''

        outstr =  BasisFn.kwrepr(self)
        
        if (outstr) or self.limitrepr():
            outstr += ','

        outstr += 'tref=' + DU.toGiantDatetime(self.tref)
        outstr += ',units=' + self.units

        return outstr
