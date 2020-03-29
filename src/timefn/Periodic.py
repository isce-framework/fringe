#!/usr/bin/env python3

############################################################
# Copyright 2013, by the California Institute of Technology#
# Author : Piyush Agram
############################################################

from .Collection import TimefnCollection
from .CenteredBasisFn import CenteredBasisFn
from .Trigonometry import Cos, Sin

class Periodic(TimefnCollection, CenteredBasisFn):
    '''
    Collection of cosine and sine of a single time period.
    '''
    
    fnname = 'Periodic'

    def __init__(self, period=None, **kwargs):
        '''
        Constructor of Periodic function of specified Time period.
        '''

        if period is None:
            raise ValueError('Time period for periodic function cannot be None')

        CenteredBasisFn.__init__(self, **kwargs)
        TimefnCollection.__init__(self)

        ###Actually build the collection
        self.build(period=period)

    def build(self, period=None):
        '''
        Build the trigonometry objects and add it to collection.
        '''

        fncos = Cos(tmin = self.tmin, tmax = self.tmax,
                    tref = self.tref, period=period,
                    units = self.units)

        fnsin = Sin(tmin = self.tmin, tmax = self.tmax,
                    tref = self.tref, period = period,
                    units = self.units)


        self.data.append(fncos)
        self.data.append(fnsin)


