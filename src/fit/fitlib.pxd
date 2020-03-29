'''
Copyright (c) 2018-, California Institute of Technology ("Caltech").
U.S. Government sponsorship acknowledged.
All rights reserved.

Author (s): Heresh Fattahi

'''

from libcpp.string cimport string

cdef extern from "Fit.h":
    cdef cppclass Fit:
        Fit() except +

        string inputDS
        string modelType;
        string outputDS;
        int memsize;
        int blocksize;
        double t_Eq;
        double Tau;
        
        void set_blocksize(int blckSize)
        void set_outputDataset(string outDS)
        void set_earthquake_time(double tEq) 
        void set_Tau(double Tau)
        
        void OpenDataset(string inputDS)
        void set_Npar(int Npar)
        void CloseDataset()
        void CreateDataset()
        void get_time()
        void fit_timeseries(string model)  

        void set_G(double* G, int n_rows, int n_cols)
        void fit_timeseries()

 
