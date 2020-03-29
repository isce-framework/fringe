#!/usr/bin/env python3

############################################################
# Copyright 2013, by the California Institute of Technology#
# Author : Piyush Agram
############################################################


from .CenteredBasisFn import CenteredBasisFn
import numpy as np 


class Trigonometry(CenteredBasisFn):
    '''Base class to represent a generic trigonometric function.

    .. math ::

        f(2 \pi \frac{t_{norm}}{T_{period}})


    Attributes
    ----------

    fnname: str
        Function name to be used for interpreting user inputs
    fnpointer: numpy function pointer
        Trigonometric numpy function to evaluate
    '''

    fnname = 'Trigonometry'

    def __init__(self, period=None, **kwargs):
        '''Constructor for Trigonometric function.

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
        period: str (or) float
            Time period for the trigonometric function
        side: str
            Alternate interface to tmin, tmax
        ''' 
        
        ####Initialize the base class
        CenteredBasisFn.__init__(self, **kwargs)
        
        try:
            self.period = float(period)
        except:
            raise TypeError('Period {0} cannot be interpreted as float'.format(period))


        if self.period == 0 :
            raise ValueError('Time period for a trigonometric function cannot be zero.')


        ###Function pointer to be set in the base class
        self.fnpointer = None
        self.radfactor = 2 * np.pi

    def computeRaw(self, tinput):
        '''Evaluation of the trigonometric function.
        '''
        return self.fnpointer( self.radfactor * self.normalizedTime(tinput)/ self.period) 

    def kwrepr(self):
        '''Keyword representation.
        '''
        outstr = CenteredBasisFn.kwrepr(self)
        if outstr:
            outstr += ','
        return  outstr + 'period=' + str(self.period)


class Cos(Trigonometry):
    '''Class to represent a cosine function.

    .. math ::

        \cos(2 \pi \frac{t_{norm}}{T_{period}})

    '''

    fnname = 'Cos'

    def __init__(self, **kwargs):
        '''Constructor for Cosine function.

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
        period: str (or) float
            Time period for the trigonometric function
        side: str
            Alternate interface to tmin, tmax
        '''

        ###Initialize the base class
        Trigonometry.__init__(self, **kwargs)

        self.fnpointer = np.cos



class Sin(Trigonometry):
    '''Class to represent a sine function.

    .. math ::

        \sin(2 \pi \frac{t_{norm}}{T_{period}})

    '''

    fnname = 'Sin'

    def __init__(self,**kwargs):
        '''Constructor for Sine function.

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
        period: str (or) float
            Time period for the trigonometric function
        side: str
            Alternate interface to tmin, tmax
        
        '''

        ###Initialize the base class
        Trigonometry.__init__(self, **kwargs)

        self.fnpointer = np.sin
   


class Arctan(Trigonometry):
    '''Class to represent an arctan function.

    .. math ::

        \tan^{-1}(\frac{t_{norm}}{T_{period}})

    '''

    fnname = 'Arctan'

    def __init__(self, **kwargs):
        '''Constructor for arctan function.

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
        period: str (or) float
            Time period for the trigonometric function
        side: str
            Alternate interface to tmin, tmax
        
        '''

        Trigonometry.__init__(self, **kwargs)
        self.fnpointer = np.arctan
        self.radfactor = 1.0
