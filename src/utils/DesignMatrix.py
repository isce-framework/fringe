#!/usr/bin/env python3

'''
Copyright (c) 2018-, California Institute of Technology ("Caltech").
U.S. Government sponsorship acknowledged.
All rights reserved.
Author (s): Heresh Fattahi
'''

import datetime
import numpy as np

class  DesignMatrix(object):
    """
    Class to create design matrix for a network of interferograms.

    Examples
    --------
    >>> import DesignMatrix
    >>> dm = DesignMatrix.DesignMatrix() 
    >>> dm.pairs=['1_2','1_3','2_3', '2_4', '2_5', '3_4', '3_5', '3_6','4_5', '4_6', '5_6']
    >>> dm.configure()
    >>> dm.closure()
    >>> G1 = dm.getG()
    >>> dm.timeseries()
    >>> G2 = dm.getG()
    >>> dm.differential_timeseries()
    >>> G3 = dm.getG()
    
    """

    def __init__(self, pairs = None, delimiter="_"):
        """Constructor for design matrix.
        
        Parameters
        ----------
        pairs: a list of strings representing reference and follower 
               acquisitions separated by delimiter which can be specified.
        delimiter: The string used to separate refernce and follower acquisitions. By default, this is '_'.

        """
        self.G = None
        self.pairs = pairs
        self.__pairs_delimiter__ = delimiter

    def configure(self, referenceEpoch = None):
        """
        configures the DesignMatrix object.
        extracts the reference and follower acquisitions and creates a list of unique epochs.

        Parameters
        ----------
        referenceEpoch: the acquisition which is considered as reference epoch for the timeseries.
        """
        self.numPairs = len(self.pairs)
        
        self.epochs = []
        self.references = [None]*self.numPairs
        self.followers = [None]*self.numPairs
        for i in range(self.numPairs):
            reference = self.pairs[i].split(self.__pairs_delimiter__)[0]
            follower = self.pairs[i].split(self.__pairs_delimiter__)[1]
            if reference not in self.epochs: self.epochs.append(reference)
            if follower not in self.epochs: self.epochs.append(follower)
            self.references[i] = reference
            self.followers[i] = follower

        self.numEpochs = len(self.epochs)
        self.epochs.sort()
        if referenceEpoch:
            self.referenceEpoch = referenceEpoch 
        else:
            self.referenceEpoch = self.epochs[0]
           
    
    def timeseries(self):
        """
        Design G matrix for a model as d = Gm where 
        d is the interferometric phases of pairs of SAR 
        acquisitions and  m is the timeseries of phase relative
        to a reference acquisition. The design matrix is same as
        A matrix in the SBAS algorithm (Berardino et al, 2002)
        
        """

        A = np.zeros((self.numPairs, self.numEpochs))
        for ni in range(self.numPairs):
            ndxt1 = self.epochs.index(self.references[ni])
            ndxt2 = self.epochs.index(self.followers[ni])
            A[ni,ndxt1] = -1
            A[ni,ndxt2] = 1
        refIndex = self.epochs.index(self.referenceEpoch)
        A = np.hstack((A[:,0:refIndex], A[:,(refIndex+1):]))
        self.G = A

    def differential_timeseries(self):
        """
        Design G matrix for a model as d = Gm where
        d is the interferometric phases of pairs of SAR
        acquisitions and  m is the differential timeseries of phase, 
        (i.e., a sequential interferometric stack). 
        The design matrix is based on (Fattahi et al, 2017).
        """
        C = np.zeros((self.numPairs, self.numEpochs-1))
        for ni in range(self.numPairs):
            ndxt1 = self.epochs.index(self.references[ni])
            ndxt2 = self.epochs.index(self.followers[ni])
            C[ni,ndxt1:ndxt2] = 1

        self.G = C

    def timeseries_timefn(self):
         self.G = None

    def closure(self):
        # copied this function from PySAR 
        """Generate the design matrix for triplets of interferograms for a model 
        where Gd = 0. 
        
        Parameters 
        ----------
        self.pairs : list of strings in t1_t2 form (t1: first acquisition, t2: second acquisition)

        Returns
        -------    
        G : 2D np.array in size of (num_tri, num_ifgram) consisting 0, 1, -1
                        for 3 SAR acquisition in t1, t2 and t3 in time order,
                        ifg1 for (t1, t2) with 1
                        ifg2 for (t1, t3) with -1
                        ifg3 for (t2, t3) with 1
        
        """
        triangle_idx = []
        for ifgram1 in self.pairs:
            # ifgram1 (t1, t2)
            t1, t2 = ifgram1.split('_')

            # ifgram2 candidates (t1, t3)
            t3_list = []
            for ifgram2 in self.pairs:
                if t1 == ifgram2.split('_')[0] and ifgram2 != ifgram1:
                    t3_list.append(ifgram2.split('_')[1])

            # ifgram2/3
            if len(t3_list) > 0:
                for t3 in t3_list:
                    ifgram3 = '{}_{}'.format(t2, t3)
                    if ifgram3 in self.pairs:
                        ifgram1 = '{}_{}'.format(t1, t2)
                        ifgram2 = '{}_{}'.format(t1, t3)
                        ifgram3 = '{}_{}'.format(t2, t3)
                        triangle_idx.append([self.pairs.index(ifgram1),
                                             self.pairs.index(ifgram2),
                                             self.pairs.index(ifgram3)])

        triangle_idx = np.array(triangle_idx, np.int16)
        triangle_idx = np.unique(triangle_idx, axis=0)

        # triangle_idx to C
        num_triangle = triangle_idx.shape[0]
        self.G = np.zeros((num_triangle, len(self.pairs)), np.float32)
        for i in range(num_triangle):
            self.G[i, triangle_idx[i, 0]] = 1
            self.G[i, triangle_idx[i, 1]] = -1
            self.G[i, triangle_idx[i, 2]] = 1

    def getG(self):
        return self.G   



