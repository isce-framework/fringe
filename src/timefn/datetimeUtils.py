#!/usr/bin/env python3

############################################################
# Copyright 2013, by the California Institute of Technology#
# Author : Piyush Agram
############################################################

import datetime
import numpy as np

DATETIME_MIN = datetime.datetime(1,1,1)  
DATETIME_MAX = datetime.datetime(3000,12,31)

####Only date
datefmt = '%Y-%m-%d'

####Date and time to 1 second
datetimefmt = '%Y-%m-%d %H:%M:%S'

####Full python date time
fulldatetimefmt = '%Y-%m-%d %H:%M:%S.%f'

def toGiantDatetime(dtobj):
    '''Convert a python datetime object into a string.

    If hrs, mins, secs and usecs are zero; only the date is returned.
    If only usecs is zero, date time without microsecs is returned.
    If not, full python datetime str is returned.

    Parameters
    ----------

    dtobj: datetime.datetime
        Python datetime object to be converted to GIAnT representation.

    '''

    if not any([dtobj.hour, dtobj.minute, dtobj.second, dtobj.microsecond]):
        return datetime.datetime.strftime(dtobj, datefmt)
    elif dtobj.microsecond == 0:
        return datetime.datetime.strftime(dtobj, datetimefmt)
    else:
        return str(dtobj)


def fromGiantDatetime(instr):
    '''Convert a string to python date time object.

    First try to parse it as a date, if not parse as datetime without microsecs.
    If this fails, try to interpret is a full python datetime object.

    Parameters
    ----------
    instr: str
        String representation of datetime in GIAnT.

    '''
    ####Try if value is in date
    try:
        result = datetime.datetime.strptime(instr, datefmt)
    except ValueError:
        result = None

    #if not result:
    if result:
        return result

    ####Try if value is in date time to second accuracy
    try:
        result = datetime.datetime.strptime(instr, datetimefmt)
    except ValueError:
        result = None

    #if not result:
    if result:
        return result

    ###Finally try, full python date time
    try:
        result = datetime.datetime.strptime(instr, fulldatetimefmt)
    except:
        raise ValueError('Could not convert {0} to GIAnT datetime type'.format(instr))

    return result
   

def interpretAsDatetime(invalue):
    '''Interpret the input correctly as a python datetime object.

    This function is meant to be used when the inputs can either be a python datetime object or a GIAnT datetime representation.

    Parameters
    ----------
    invalue: datetime.datetime (or) datetime.date (or) str
        Input to be interpreted as python datetime

    Returns
    -------
    datetime.datetime
        Input translated into a python datetime object
    '''
    if isinstance(invalue, str):
        outvalue = fromGiantDatetime(invalue) # tmin)
    elif isinstance(invalue, datetime.datetime):
        outvalue = invalue
    elif isinstance(invalue, datetime.date):
        outvalue = datetime.datetime.fromordinal(invalue.toordinal())
    else:
        raise ValueError('Datetime values can only be provided as GIAnT compatible strings or python datetime objects')

    return outvalue


def secondsOfDay(invalue):
    '''Returns seconds of day as float.

    Parameters
    ----------
    invalue : str (or) datetime.date (or) datetime.datetime
        Any datetime input handled by GIAnT

    Returns
    -------
    float
        Seconds of day as a float

    '''

    dt = interpretAsDatetime(invalue)
    midnight = datetime.datetime.combine(dt.date(), datetime.time(0))

    return (dt - midnight).total_seconds()


def hoursOfDay(invalue):
    '''Returns hours of day as float.

    Parameters
    ----------
    invalue : str (or) datetime.date (or) datetime.datetime
        Any datetime input handled by GIAnT

    Returns
    -------
    float
        Hours of day as a float
    '''

    secs = secondsOfDay(invalue)
    return secs / 60.0

def dayOfYear(invalue):
    '''Returns day of year as float.

    Parameters
    ----------
    invalue : str (or) datetime.date (or) datetime.datetime
        Any datetime input handled by GIAnT

    Returns
    -------
    float
        Day of year as a float
    '''

    dt = interpretAsDatetime(invalue)
    int_day = dt.timetuple().tm_yday

    return int_day + secondsOfDay(dt) / 86400.0


def linspace(tmin, tmax, num=None):
    '''
    Equivalent to numpy linspace but for python datetime objects.

    Parameters
    ----------
    tmin: str (or) datetime.date (or) datetime.datetime
        Any GIAnT compatible datetime format
    tmax: str (or) datetime.date (or) datetime.datetime
        Any GIAnT compatible datetime format
    num: int
        Number of samples to divide the interval into


    Returns
    -------
    np.array
        Array of python datetime objects

    '''

    outlist = []

    tmin = interpretAsDatetime(tmin)
    tmax = interpretAsDatetime(tmax)

    try:
        num = int(num)
    except:
        raise ValueError('datetime linspace needs integer number of samples. {0} cannot be interpreted as an int'.format(num))

    dt = (tmax - tmin).total_seconds()/(num-1.0)

    for x in range(num):
        t = tmin + datetime.timedelta(seconds = x * dt)
        outlist.append(t)

    return np.array(outlist)
