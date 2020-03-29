#ifndef  Fit_timeseries_h
#define  Fit_timeseries_h

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

class Fit {
  public:

    Fit(void) {};
    ~Fit(void) {};

    static Fit* instance(void)
            {
              if (self == NULL) {self = new Fit();}
              return self;
            };
    

    virtual void OpenDataset(str inputDS);
    virtual void OpenDataset(str inputDS, str tempCorDSName);
    virtual void CloseDataset();

    virtual void CreateDataset();
    virtual void CreateDataset(str outputDS);
    virtual void CreateDataset(str outputDS, int xSize, int ySize, int zSize);
    virtual void CreateDataset(str outputDS, int xSize, int ySize , int zSzie, str driverName);

    virtual void get_shape(GDALDataset* Dataset, int &cols, int &rows, int &nbands);

    //virtual void define_blocksize();

    virtual void get_time();
    //virtual void get_data_block(int ysize, int yoff, arma::mat& data);
    virtual void get_data_block(int ysize, int yoff, arma::Mat<float>& data);
    virtual void get_data_block(int ysize, int yoff, arma::Mat<float>& data, arma::Mat<float>& tempCoh);
    virtual void define_design_matrix(arma::mat& G, arma::mat& G_inv, str model);


    virtual void set_blocksize(int blckSize);
    virtual void set_Npar(int par);
    virtual void set_outputDataset(str outDS);
    virtual void set_tauDataset(str tauDS);
    virtual void set_earthquake_time(double tEq);    
    virtual void set_Tau(double tau);
    virtual void set_temporal_coherence_threshold(float tcor);

    virtual void Heaviside(arma::mat& H);
    virtual void PostSeismic(arma::mat& P, arma::mat& H, double Tau);

    virtual void set_G(double* G_temp, int n_rows, int n_cols);
    virtual void fit_timeseries();

    virtual void fit_timeseries(str model);
    virtual void fit_timeseries(str model, float &meanRMSE);
    virtual void write2datset(arma::Mat<float> X, int Npar, int yoff, int ysize);
    virtual void write2datset(arma::mat X, arma::mat residual, int Npar, int yoff, int ysize);   

    virtual void estimate_Tau(str model);
    virtual void get_sum(arma::mat rmse, arma::Mat<float> tcorr, float threshold , float &sum, int &count);
    virtual void write_Tau(arma::mat X);
    //virtual void linear_design_matrix(arma::mat& B);
    //virtual void set_design_matrix(str model, float Tau);
    //virtual void solve();
    //virtual void get_residual();
    //virtual void write_parameters();
 
    str inputDS;
    str modelType;
    str outputDS;
    str tauDS;
    str tempCorrDS;
    int memsize;
    int blocksize;
    int rows, cols, nbands; // size of input dataset
    int Npar; // number of unknown parameters
    double t_Eq; // the time of an earthquake
    double Tau; // post-seismic relaxation time
    float tcorrThreshold; // threshold of temporal coherence
    arma::mat G;

    private:
        static Fit* self;

        Fit(const Fit&) {assert(false);}

        void operator=(const Fit&) {assert(false);}

        // private data members
        GDALDataset* inDataset = NULL;
        GDALDataset* outDataset = NULL;
        GDALDataset* tcorDataset = NULL;
        GDALDataset* tauDataset = NULL;
        arma::mat acquisition_times;
        /*arma::mat reference_phase; //(nbands, 1);
        arma::mat Heaviside; // (nbands, 1);
        arma::mat A; //  (nbands, 3);
        arma::mat A_inv;
        arma::mat X;
        arma::mat rates;
        arma::mat coseismic;*/

}; 
     
#endif

