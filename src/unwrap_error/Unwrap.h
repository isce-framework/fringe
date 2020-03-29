#ifndef  Unwrap_timeseries_h
#define  Unwrap_timeseries_h

#include "gdal.h"
#include "gdal_priv.h"
#include <vector>
#include <armadillo>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <iomanip>
#include <cassert>

double SPEED_OF_LIGHT = 299792458.0;
const double PI = 3.141592653589793;


typedef std::string str;

using namespace std;

class Unwrap {
   
   public:

    Unwrap(void) {};
    ~Unwrap(void) {};

    static Unwrap* instance(void)
            {
              if (self == NULL) {self = new Unwrap();}
              return self;
            };

    virtual void Set_connComp_dataset(str inDS);
    virtual void Set_mask_file(str maskFile);

    virtual void OpenDataset();
    virtual void CreateMaskDataset();
    virtual void CreateMaskDataset(str outputDS, int xSize, int ySize , int zSize, str driverName, GDALDataType dtype);
    virtual void CloseDataset();

    virtual void ComputeCommonMask();
     
    str connCompDS; // connected components data set
    str maskDS;     // A mask dataset obtained from multiplication of all conn components
    int blocksize;
    int rows, cols, nbands; // size of input dataset
    
    private:
        static Unwrap* self;

        Unwrap(const Unwrap&) {assert(false);}

        void operator=(const Unwrap&) {assert(false);}

        GDALDataset* inDataset = NULL;
        GDALDataset* maskDataset = NULL;
        virtual void compute_mask(int, int, arma::Mat<short int>&);
        virtual void write2datset(arma::Mat<short>, int, int );

};

#endif
