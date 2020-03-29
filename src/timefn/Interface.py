#!/usr/bin/env python3

############################################################
# Copyright 2013, by the California Institute of Technology#
# Author : Piyush Agram
############################################################

####Collection of functions in one place
from collections import OrderedDict

def addFn(rdict, fnlist):

    ###For each function in input list
    for fn in fnlist:

        ###Get name of function which is class variable
        key = fn.fnname

        ####Check if key already exists in function mapping
        if key in rdict.keys():
            raise Exception('Function collection already contains a function with key {0}'.format(key))

        ####Normalize key using lower
        rdict[key.lower()] = fn
    return

fnmap = OrderedDict()   ###Mapping of fnnames to classes
regmap = OrderedDict() ###Mapping of fnnames to regularization booleans

###Step function
from .Constant import Constant
addFn(fnmap, [Constant])

####Power functions
from .Power import Power, Linear, Quadratic, Cubic, Quartic
addFn(fnmap, [Power, Linear, Quadratic, Cubic, Quartic])

####Log functions
from .Log import Log, Log10, Log2, GeoLog
addFn(fnmap, [Log, Log10, Log2, GeoLog])

####Exp functions
from .Exp import Exp, Exp10, Exp2
addFn(fnmap, [Exp, Exp10, Exp2])

####Trigonometry functions
from .Trigonometry import Cos, Sin, Arctan
addFn(fnmap, [Cos, Sin, Arctan])

####Spline functions
from .Spline import BSpline, ISpline, PBSpline
addFn(fnmap, [BSpline, ISpline, PBSpline])

####Polynomial
from .Polynomial import Polynomial
addFn(fnmap, [Polynomial])

####Periodic
from .Periodic import Periodic
addFn(fnmap, [Periodic])

###Spline sets
from .SplineFactory import BSplineSet, ISplineSet, PBSplineSet
addFn(fnmap, [BSplineSet, ISplineSet, PBSplineSet])


###Collection
from .Collection import TimefnCollection


def fromString(instr):
    '''
    This function is responsible for taking a string and interpreting it as a collection of functions (or) individual function.

    Parameters
    ----------
    instr: str
        string to be interpreted as a timefn function.

    Returns
    -------
    obj
        python timefn/ timefncollection object corresonding to the string

    '''
    
    ####Kicks out spaces from beginning / end of string
    instr = instr.strip()

    
    if instr.startswith('['):
        ####Translate string to collection
        result = collectionFromString(instr)
    else:
        ####Translate string to function
        result = timefnFromString(instr)

    return result


def timefnFromString(instr):
    '''
    This function is responsible for taking in a string and interpreting as an individual function.

    Parameters
    ----------
    instr: str
        string to be interpreted as a timefn function.

    Returns
    -------
    obj
        python timefn object corresponding to the string.

    '''

    ###Utility functions
    def genError(instring):
        raise ValueError('Could not interpret string "{0}" as timefn'.format(instring))

    def genPairError(instring):
        raise ValueError('Could not interpret string "{0}" as kwarg'.format(instring))


    ###Kicks out spaces from beginning/ end of string
    instr = instr.strip()
    ###Split with '(' to get the function name and the rest
    try:
        fnname, rest = instr.split('(')#[0]
    except:
        genError(instr)
   
    ###Parse the rest of the string.
    if not rest.endswith(')'):
        genError(instr)

    ####Kick out the ')' at the end
    rest = rest[:-1]

    ####Split what is left into kwarg pairs
    fields = rest.split(',')

    ###Empty dictionary
    kwdict = {}

    ###For each pair
    for field in fields:
        try:
            key, val = field.split('=')
        except:
            genPairError(field)
       
        ###Take care of any quotes that the user may have put in
        val = val.replace("'",'')
        val = val.replace('"','')

        kwdict[key.strip()] = val.strip()

    ###Lookup function name in fnmap
    try:
        basefn = fnmap[fnname.lower()]
    except:
        raise KeyError('Could not find function named {0}'.format(fnname.lower()))


    return basefn(**kwdict)


def collectionFromString(instr):
    '''
    This function is responsible for taking in a string and interpreting it as a collection of timefn.

    Parameters
    ----------
    instr: str
        string to be interpreted as a timefn function.

    Returns
    -------
    obj
        python timefncollection object corresponding to the string.

    '''
    
    ###Utility functions
    def genError(instring):
        raise ValueError('Cound not interpret string "{0}" as a timefn collection. Collections should be contained with [ ] and functions must be separated by ;.'.format(instring))

    ###Kick out empty spaces at the start/end of string
    instr = instr.strip()

    ####Create an empty collction
    collection = TimefnCollection()

    ###See if string is contained within square brackets
    if instr.startswith('[') and instr.endswith(']'):
        fns = instr[1:-1].split(';')
        
        for fn in fns:
            newfn = timefnFromString(fn)
            collection.append(newfn)
            #collection += newfn
    else:
        genError(instring)

    return collection


def getFunctionTypes(collection):
    """
    Iterate over a collection to get the function types.

    Parameters
    ----------
    collection: TimefnCollection
        Input TimefnCollection object.

    Returns
    -------
    dict
        Dictionary containing indices for each function type.
    """
    secular = []; seasonal = []; transient = []; step = []
    reg = []
    for cnt, basis in enumerate(collection.data):

        if basis.fnname in ['Constant', 'Linear']:
            secular.append(cnt)
        elif basis.fnname in ['Cos', 'Sin']:
            seasonal.append(cnt)
        elif basis.fnname in ['Step']:
            step.append(cnt)
        elif basis.fnname in ['BSpline']:
            seasonal.append(cnt)
            reg.append(cnt)
        elif basis.fnname in ['ISpline']:
            transient.append(cnt)
            reg.append(cnt)

    return {'secular': secular, 'seasonal': seasonal, 'transient': transient,
            'step': step, 'reg': reg}


def getTimescales(t, collection):
    """
    Return a list of the effective timescales of each basis function
    in a collection given an array of times t.

    Parameters
    ----------
    t: list, ndarray
        Array of datetime.date/datetime.datetime objects. 
    collection: TimefnCollection
        Input TimefnCollection object

    Returns
    -------
    ndarray
        Array of timescales for each basis function in collection.
    """
    import numpy as np
    scales = []
    t0, tf = t[0], t[-1]
    tspan = abs((tf - t0).days) / 365.0
    for cnt, basis in enumerate(collection.data):
        if basis.fnname in ['Constant', 'Linear']:
            scales.append(tspan)
        elif basis.fnname in ['Step']:
            scales.append(1.0/365.0)
        elif basis.fnname in ['Cos', 'Sin']:
            scales.append(basis.period)
        elif basis.fnname in ['BSpline', 'ISpline']:
            scales.append(basis.scale)
        else:
            scales.append(np.nan)

    return np.array(scales) 

# end of file
