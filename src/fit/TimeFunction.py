#!/usr/bin/env python3

'''
Copyright (c) 2018-, California Institute of Technology ("Caltech").
U.S. Government sponsorship acknowledged.
All rights reserved.
Author (s): Heresh Fattahi
'''

import datetime
from timefn import fnmap, linspace, fromString
import numpy as np

class TimeFunction(object):

    def __init__(self, G=None, time_array = None,  units='days'):
        self.G = G
        self.time_array = time_array
        self.units = units
        self.timeFnString = ""

    def __updateFnString__(self, st):

        self.timeFnString+= st + " + "

    def configure(self):
        self.epoch_numbers = len(self.time_array)
        self.tmin = self.time_array[0]
        self.tmax = self.time_array[-1]
        if self.G is None:
            self.G = np.ones((self.epoch_numbers,1))
            self.__updateFnString__("c")
        
    def create_time_array(self, tmin , tmax, number_of_epochs):
        self.time_array = linspace(tmin, tmax, num=number_of_epochs)

    def add_heaviside(self, tmin):
        objConst = fromString("[Constant(tmin={0})]".format(tmin))
        G_heaviside = objConst(self.time_array)
        self.G = np.hstack((self.G, G_heaviside.reshape(self.epoch_numbers,1)))
        self.__updateFnString__("H({0}) ".format(tmin))

    def add_periodic(self,tref, period):
        base = fnmap['periodic']
        objPeriodic = base(tref=tref, units=self.units, period=period,
                           tmin=self.tmin, tmax=self.tmax)

        Gperiodic = objPeriodic(self.time_array)
        self.G = np.hstack((self.G, Gperiodic))
        self.__updateFnString__("sin(2*pi*t/{0}) + cos(2*pi*t/{0})".format(period))

    def add_post_seismic(self, tref, tau):
        base = fnmap['geolog']
        obj = base(tref=tref, units=self.units, tau=tau)
        Glog = obj(self.time_array)
        Glog = Glog.reshape(self.epoch_numbers,1)
 
        objConst = fromString("[Constant(tmin={0})]".format(tref))
        G_heaviside = objConst(self.time_array)

        Glog = G_heaviside*Glog

        self.G = np.hstack((self.G, Glog.reshape(self.epoch_numbers,1)))   
        self.__updateFnString__("H*(log(1.0+(t-{0})/{1}))".format(tref, tau))

    def add_log(self, tref, tau):
        base = fnmap['log']
        obj = base(tref=tref, units=self.units, tau=tau)
        Glog = obj(self.time_array)
        
        self.G = np.hstack((self.G, Glog.reshape(self.epoch_numbers,1))) 
        self.__updateFnString__("H*(log((t-{0})/{1}))".format(tref, tau))

    def add_log10(self, tref, tau):
        base = fnmap['log10']
        obj = base(tref=tref, units=self.units, tau=tau)
        Glog = obj(self.time_array)
        self.G = np.hstack((self.G, Glog.reshape(self.epoch_numbers,1)))
        self.__updateFnString__("H*(log10((t-{0})/{1}))".format(tref, tau))

    def add_polynomial(self, tref, order, tau=1):
        base = fnmap['poly']
        obj = base(tref=tref, units=self.units, tau=tau, order=order)
        Gpoly = obj(self.time_array)
        self.G = np.hstack((self.G, Gpoly[:,1:]))
        for i in range(order):
            self.__updateFnString__("(t-{0})^{1}".format(tref,i+1))
 
    def getG(self):
        return self.G, self.timeFnString[:-2]

