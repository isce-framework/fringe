#!/usr/bin/env python3

############################################################
# Copyright 2018, by the California Institute of Technology#
# Author : Piyush Agram
############################################################


from . import datetimeUtils as DU
from datetime import datetime
import numpy as np

class BasisFn(object):
    '''Base class for timefn representation.

    All timefn types are derived from this elementary type. This function only implements the active time period - minimum and maximum time limits for a function to be active.

    Attributes
    ----------

    tmin : datetime.datetime
        Minimum time
    tmax : datetime.datetime
        Maximum time
    fnname : str
        Unique name for each function.

    '''

    fnname = 'Basis'  #Unique name for function

    def __init__(self, tmin=DU.DATETIME_MIN, tmax=DU.DATETIME_MAX,  **kwargs):
        '''Constructor of the BasisFn class.

        Parameters
        ----------
        tmin : str (or) datetime.datetime
            Left edge of active period of function (keyword)
        tmax : str (or) datetime.datetime
            Right edge of active period of function (keyword)

        '''

        self.tmin = DU.interpretAsDatetime(tmin)

        self.tmax = DU.interpretAsDatetime(tmax)

        ###Use numpy vectorize to handle scalar or array
        self.isActive = np.vectorize(self.isActiveOne)

    
    def isActiveOne(self, tinput):
        '''Checks if function is active on a given input date.

        Parameters
        ----------
        tinput : datetime.datetime
            Datetime object to see if function is active

        Returns
        -------
        bool
            True if active, else False
        '''

        if (self.tmin <= tinput) and (self.tmax >= tinput):
            return 1
        else:
            return 0

    
    def computeRaw(self, tinput):
        '''Compute the functional form without applying isActive mask.

        Parameters
        ----------
        tinput: np.array of datetime.datetime
            Array of datetime objects for which the function needs to be evaluated.

        Returns
        -------
        array: np.array of type double
            Array of function values at datetimes corresponding to tinput.

        '''

        raise NotImplementedError('Specific implementation in derived classes')

    def compute(self, tinput):
        '''Computes value and applies the isActive mask.

        Parameters
        ----------
        tinput: np.array of datetime.datetime
            Array of datetime objects for which the function needs to be evaluated.

        Returns
        -------
        array : np.array of type double
            Array of function values at datetimes corresponding to tinput.
        '''

        return self.computeRaw(tinput) * self.isActive(tinput)

    def __call__(self, tinput):
        ''' Overloading of the compute function.

        Note
        -----
        Overloading the "compute" method. Meant to be the user interface in the code.

        '''
        return self.compute(tinput)


    def __repr__(self):
        '''String representation of basis function.

        Note
        ----
        Returns the string representation of a basis function.
        The string representation consists of 3 parts:
        *   Function name followed by '('
        *   Keyword and value pairs, separated by a comma.
        *   Terminated by ')'

        '''

        return self.preamble() + self.limitrepr() + self.kwrepr() + ')'

    def __eq__(self, other):
        '''Compare for equality.

        Note
        ----

        Simple comparison is performed by checking string representations.
        '''
        
        return repr(other) == repr(self)

    def preamble(self):
        '''String representation of function name.

        Note
        ----

        Implments the first part of the string representation.
        '''
        return self.fnname+'('

    def limitrepr(self):
        '''String representation of the active period.

        Notes
        -----

        Sets up the active period in the string representation.
        '''
        outstr = ''

        if (self.tmin != DU.DATETIME_MIN):
            outstr += 'tmin=' + DU.toGiantDatetime(self.tmin)

        if (self.tmax != DU.DATETIME_MAX):
            if outstr:
                outstr += ','
            outstr += 'tmax=' + DU.toGiantDatetime(self.tmax)

        return outstr

    def kwrepr(self):
        '''Specific for representing keywords.

        Note
        ----

        This function will be customized for each implemented function.
        '''
        return ''




