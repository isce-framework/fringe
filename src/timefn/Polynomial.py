#!/usr/bin/env python3

############################################################
# Copyright 2013, by the California Institute of Technology#
# Author : Piyush Agram
############################################################

from .Collection import TimefnCollection
from .CenteredBasisFn import CenteredBasisFn
from .Power import Power
from .Constant import Constant
from .Power import Linear, Quadratic, Cubic, Quartic

fmap = { 1 : Linear,
         2 : Quadratic,
         3 : Cubic,
         4 : Quartic}

class Polynomial(TimefnCollection, CenteredBasisFn):
    '''
    Collection of power objects to represent a polynomial.
    '''

    fnname = 'Poly'

    def __init__(self, order=None, tau=1, minorder=0, **kwargs):
        '''
        Constructor of a polynomial of specified order.
        '''

        if order is None:
            raise ValueError('Polynomial order cannot be None')

        CenteredBasisFn.__init__(self, **kwargs)
        TimefnCollection.__init__(self)

        ###Actually build the collection
        self.build(order=order, tau=tau, minorder=minorder)

    def build(self, order=None, tau=1, minorder=0):
        '''
        Build the power objects and add it to collection.
        '''

        for exp in range(minorder, order+1):
           
            if exp == 0:
                fn = Constant(tmin=self.tmin, tmax=self.tmax)

            elif (exp in fmap.keys()):
                fn  = fmap[exp](tmin = self.tmin, tmax = self.tmax,
                                tref = self.tref, tau = tau,
                                units = self.units)
            else:
                fn  = Power(tmin = self.tmin, tmax = self.tmax,
                            tref = self.tref, tau = tau,
                            units = self.units, exp=exp)

            self.data.append(fn)


