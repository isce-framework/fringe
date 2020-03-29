#!/usr/bin/env python3

############################################################
# Copyright 2013, by the California Institute of Technology#
# Author : Piyush Agram
############################################################

from .Collection import TimefnCollection
from .BasisFn import BasisFn
from .Spline import BSpline, ISpline, PBSpline #, PBSpline
from .datetimeUtils import linspace
from .CenteredBasisFn import ureg, CenteredBasisFn 
from . import datetimeUtils as DU
import numpy as np 

class SplineFactory(TimefnCollection, BasisFn):
    '''Base class to construct a collection of spline objects.

    Attributes
    ----------
    tmin : str (or) datetime.date (or) datetime.datetime
        Any Giantdatetime input. Left edge of active period.
    tmax : str (or) datetime.date (or) datetime.datetime
        Any Giantdatetime input. Right edge of active period.
    order : int
        Order for splines to build
    num : int
        Number of splines to build
    '''

    fnname = 'SplineSet'

    def __init__(self, order=2, num=None, units='secs', **kwargs):
        '''Constructor of a uniform set of splines

        Parameters
        ----------
        tmin : str (or) datetime.date (or) datetime.datetime
            Any Giantdatetime input. Left edge of active period.
        tmax : str (or) datetime.date (or) datetime.datetime
            Any Giantdatetime input. Right edge of active period.
        order : int
            Order for splines to build
        num : int
            Number of splines to build
        units : str
            Python pint compatible units string
        '''

        try: 
            order = int(order)
        except:
            raise ValueError('Spline order cannot be interpreted as an int')

        if order == 0 :
            raise ValueError('Spline order cannot be zero.')

        try:
            num = int(num)
        except:
            raise ValueError('Number of splines cannot be interpreted as int')

        if (num == 0):
            raise ValueError('Number of splines cannot be zero')

        BasisFn.__init__(self, **kwargs)
        TimefnCollection.__init__(self)
        self.fnname = 'SplineFactory'

        #self.baseClass = None

        ###Actually build the collection
        self.build(order=order, num=num, units=units)

    def build(self, order=None, num=None, units=None):
        '''Build the spline objects and add it to collection.

        Uses the baseClass pointer to build a spline of a particular kind.
        '''

        if units is None:
            raise ValueError('No units provided for building splines') 

        knots = DU.linspace(self.tmin, self.tmax, num=num)
        scale = (knots[1] - knots[0]).total_seconds()

        ####Convert scale to same units as requested by user
        scale_units = (scale*ureg['secs']).to(units).magnitude

        for knot in knots:
            fn = self.baseClass(tref=knot, order=order,
                                units=units, scale=scale_units)

            self.data.append(fn)


class BSplineSet(SplineFactory):
    '''A collection of bspline objects.

    Attributes
    ----------
    tmin : str (or) datetime.date (or) datetime.datetime
        Any Giantdatetime input. Left edge of active period.
    tmax : str (or) datetime.date (or) datetime.datetime
        Any Giantdatetime input. Right edge of active period.
    order : int
        Order for splines to build
    num : int
        Number of splines to build
    '''

    fnname = 'BSplineSet'

    def __init__(self, **kwargs):
        '''Constructor of a uniform set of bsplines

        Parameters
        ----------
        tmin : str (or) datetime.date (or) datetime.datetime
            Any Giantdatetime input. Left edge of active period.
        tmax : str (or) datetime.date (or) datetime.datetime
            Any Giantdatetime input. Right edge of active period.
        order : int
            Order for splines to build
        num : int
            Number of splines to build
        units : str
            Python pint compatible units string
        '''
        self.baseClass = BSpline
        SplineFactory.__init__(self, **kwargs)


class ISplineSet(SplineFactory):
    '''A collection of ispline objects.

    Attributes
    ----------
    tmin : str (or) datetime.date (or) datetime.datetime
        Any Giantdatetime input. Left edge of active period.
    tmax : str (or) datetime.date (or) datetime.datetime
        Any Giantdatetime input. Right edge of active period.
    order : int
        Order for splines to build
    num : int
        Number of splines to build
    '''

    fnname = 'ISplineSet'

    def __init__(self, **kwargs):
        '''Constructor of a uniform set of isplines

        Parameters
        ----------
        tmin : str (or) datetime.date (or) datetime.datetime
            Any Giantdatetime input. Left edge of active period.
        tmax : str (or) datetime.date (or) datetime.datetime
            Any Giantdatetime input. Right edge of active period.
        order : int
            Order for splines to build
        num : int
            Number of splines to build
        units : str
            Python pint compatible units string
        '''        
        self.baseClass = ISpline
        SplineFactory.__init__(self, **kwargs)
        

class PBSplineSet(TimefnCollection, CenteredBasisFn):
    '''A collection of pbspline objects.

    Attributes
    ----------
    tmin : str (or) datetime.date (or) datetime.datetime
        Any Giantdatetime input. Left edge of active period.
    tmax : str (or) datetime.date (or) datetime.datetime
        Any Giantdatetime input. Right edge of active period.
    order : int
        Order for splines to build
    num : int
        Number of splines to build
    period : float
        Time period for the splines
    '''

    fnname = 'PBSplineSet'

    def __init__(self, order=2, num=None,
                       period=None,  **kwargs):
        '''Constructor of a uniform set of splines

        Parameters
        ----------
        tmin : str (or) datetime.date (or) datetime.datetime
            Any Giantdatetime input. Left edge of active period.
        tmax : str (or) datetime.date (or) datetime.datetime
            Any Giantdatetime input. Right edge of active period.
        order : int
            Order for splines to build
        num : int
            Number of splines to build
        period : float
            Time period for the splines
        units : str
            Python pint compatible units string
        '''

        try: 
            order = int(order)
        except:
            raise ValueError('Spline order {0} cannot be interpreted as an int'.format(order))

        if order == 0 :
            raise ValueError('Spline order cannot be zero.')

        try:
            num = int(num)
        except:
            raise ValueError('Number of splines {0} cannot be interpreted as int'.format(num))

        if (num == 0):
            raise ValueError('Number of splines cannot be zero')

        try:
            period = float(period)
        except:
            raise ValueError('Time period {0} for the splines cannot be interpreted as float'.format(period))

        if (period == 0):
            raise ValueError('Time period for the splines cannot be zero')

        CenteredBasisFn.__init__(self, **kwargs)
        TimefnCollection.__init__(self)

        self.baseClass = PBSpline

        ###Actually build the collection
        self.build(order=order, num=num, units=units, period=period)


    def build(self, order=None, num=None, units=None, period=None):
        '''Build the spline objects and add it to collection.

        Uses the baseClass pointer to build a spline of a particular kind.
        '''

        if units is None:
            raise ValueError('No units provided for building splines') 

        knots = np.linspace(0., period, num=num+1)
        scale = knots[1] - knots[0]

        for knot in knots[:-1]:
            fn = self.baseClass(tref=self.tref, order=order,
                                units=units, scale=scale,
                                offset=knot, period=period)

            self.data.append(fn)


