'''
Copyright (c) 2018-, California Institute of Technology ("Caltech").
U.S. Government sponsorship acknowledged.
All rights reserved.

Author (s): Heresh Fattahi

'''

from libcpp.string cimport string

cdef extern from "Unwrap.h":
    cdef cppclass Unwrap:
        Unwrap() except +

        string connCompDS
        string maskDS
        int blocksize
        
        void Set_connComp_dataset(string inputDS)
        void Set_mask_file(string outDS)
        
        void OpenDataset()
        void CloseDataset()
        void CreateMaskDataset()


        void ComputeCommonMask()

        #void set_G(double* G, int n_rows, int n_cols)
        #void estimate_unwrap_error() 
         
 
