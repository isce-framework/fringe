// ----------------------------------------------------------------------------
// Copyright (c) 2017-, California Institute of Technology ("Caltech"). U.S.
// Government sponsorship acknowledged.
// All rights reserved.
//  
// Author(s): Heresh Fattahi
// ----------------------------------------------------------------------------

#include "Fit.h"



void Fit::OpenDataset(str inputDS){

  GDALAllRegister();
  this->inDataset = reinterpret_cast<GDALDataset *>( GDALOpenShared( inputDS.c_str(), GA_ReadOnly));    

  if (inDataset == NULL)
    {
        std::cout << "Cannot open stack file { " <<  inputDS << " } with GDAL for reading. \n";

        std::cout << "GDALOpen failed - " << inputDS << "\n";
        std::cout << "Exiting with error code .... (102) \n";
        GDALDestroyDriverManager();
        //return 102;
    }

    this->cols = inDataset->GetRasterXSize();
    this->rows = inDataset->GetRasterYSize();
    this->nbands = inDataset->GetRasterCount();


}


void Fit::OpenDataset(str inputDS, str tempCorDSName){

  GDALAllRegister();
  this->inDataset = reinterpret_cast<GDALDataset *>( GDALOpenShared( inputDS.c_str(), GA_ReadOnly));
  
  this->cols = inDataset->GetRasterXSize();
  this->rows = inDataset->GetRasterYSize();
  this->nbands = inDataset->GetRasterCount();

  this->tcorDataset = reinterpret_cast<GDALDataset *>( GDALOpenShared( tempCorDSName.c_str(), GA_ReadOnly));

}

void Fit::CloseDataset(){

    std::cout << "closing the dataset" << std::endl;
    if (inDataset != NULL)
        GDALClose(inDataset);

    if (outDataset != NULL)
        GDALClose(outDataset);

    if (tcorDataset !=NULL)
        GDALClose(tcorDataset);

    std::cout << "Destroy the driver" << std::endl;
    GDALDestroyDriverManager();
    

}

void Fit::set_outputDataset(str outDS){

   this->outputDS=outDS;

}

void Fit::set_tauDataset(str tauDS){

   this->tauDS=tauDS;

}

void Fit::CreateDataset(){

    str driverName = "ENVI";
    CreateDataset(outputDS,  cols, rows , Npar, driverName);

}


void Fit::CreateDataset(str outDS){

    this->outputDS=outDS;
    str driverName = "ENVI";
    CreateDataset(outputDS,  cols, rows , Npar, driverName);

}

void Fit::CreateDataset(str outputDS, int xSize, int ySize , int zSize){

    str driverName = "ENVI";
    CreateDataset(outputDS,  xSize, ySize , zSize, driverName);

}

void Fit::CreateDataset(str outputDS, int xSize, int ySize , int zSize, str driverName){

    GDALAllRegister();
    GDALDriver *poDriver = (GDALDriver*) GDALGetDriverByName(driverName.c_str());
    char **mOptions = NULL;
    mOptions = CSLSetNameValue(mOptions, "INTERLEAVE", "BIP");
    mOptions = CSLSetNameValue(mOptions, "SUFFIX", "ADD");
    outDataset = (GDALDataset*) poDriver->Create(outputDS.c_str(), xSize, ySize, zSize, 
                                                GDT_Float32, mOptions);
    if (outDataset == NULL)
    {
        std::cout << "Could not create output dataset {" << outputDS << "} \n";
        std::cout << "Exiting with non-zero error code ... 104 \n";

        GDALClose(inDataset);
        GDALClose(outDataset);
        GDALDestroyDriverManager();
    }

    CSLDestroy(mOptions);

}



void Fit::write2datset(arma::Mat<float> X, int Npar, int yoff, int ysize){
    for (int ii =0; ii<Npar; ii++){
        int status = outDataset->GetRasterBand(ii+1)->RasterIO(GF_Write, 0, yoff,
                    cols, ysize,
                    (void*) (&X.col(ii)[0]),
                    cols, ysize, GDT_Float32,
                    0 , 0, NULL);

                   // sizeof(double), sizeof(double)*cols, NULL);
    
    }
}

void Fit::write2datset(arma::mat X, arma::mat residual, int Npar, int yoff, int ysize){
    for (int ii =0; ii<Npar; ii++){
        int status = outDataset->GetRasterBand(ii+1)->RasterIO(GF_Write, 0, yoff,
                    cols, ysize,
                    (void*) (&X.col(ii)[0]),
                    cols, ysize, GDT_Float32,
                    0 , 0, NULL);

    }

    int status = outDataset->GetRasterBand(Npar+1)->RasterIO(GF_Write, 0, yoff,
                    cols, ysize,
                    (void*) (&residual[0]),
                    cols, ysize, GDT_Float32,
                    0 , 0, NULL);

}

 
void Fit::get_shape(GDALDataset *Dataset, int &cols, int &rows, int &nbands){ 
  cols = Dataset->GetRasterXSize();
  rows = Dataset->GetRasterYSize();
  nbands = Dataset->GetRasterCount();

  std::cout << "Number of rows  = " << rows << "\n";
  std::cout << "Number of cols  = " << cols << "\n";
  std::cout << "Number of bands = " << nbands << "\n";
  
}


void Fit::set_blocksize(int blckSize){

    this-> blocksize = blckSize;

}


/*void Fit::define_blocksize(){
  blocksize = int((memsize * 1.0e6)/(cols * 8 * 6) );
  if (blocksize > rows)
    blocksize = rows;
  std::cout << "Computed block size based on memory size = " << blocksize << " lines \n";
  blocksize = 128;

}*/

void Fit::get_time(){

   acquisition_times.set_size(nbands,1);
   double acquisition;
  for (int bb=0; bb< nbands; bb++){
     acquisition = std::stod(inDataset->GetRasterBand(bb+1)->GetMetadataItem("AcquisitionTime", "slc"));
     acquisition_times(bb,0) = acquisition;
     

     }
    
    //acquisition_times.print();


}


//void Fit::get_data_block(int ysize, int yoff, arma::mat &data){
void Fit::get_data_block(int ysize, int yoff, arma::Mat<float> &data){
   data.set_size(cols*ysize, nbands);

   for (int bb=0; bb< nbands; bb++){
            int status = inDataset->GetRasterBand(bb+1)->RasterIO( GF_Read, 0, yoff,
                                cols, ysize,
                                (void*) (data.colptr(bb)),
                                cols, ysize, GDT_Float32,
                                sizeof(float),
                                sizeof(float)*cols, NULL);
            //data.col(bb) = data.col(bb); //+ reference_phase(bb,0);
        }

}

void Fit::get_data_block(int ysize, int yoff, arma::Mat<float> &data, arma::Mat<float> &tempCoh){

    data.set_size(cols*ysize, nbands);
    tempCoh.set_size(cols*ysize, 1);   

    for (int bb=0; bb< nbands; bb++){
            int status = inDataset->GetRasterBand(bb+1)->RasterIO( GF_Read, 0, yoff,
                                cols, ysize,
                                (void*) (data.colptr(bb)),
                                cols, ysize, GDT_Float32,
                                sizeof(float),
                                sizeof(float)*cols, NULL);
        }

    int status = tcorDataset->GetRasterBand(1)->RasterIO( GF_Read, 0, yoff,
                                cols, ysize,
                                (void*) (tempCoh.colptr(0)),
                                cols, ysize, GDT_Float32,
                                sizeof(float),
                                sizeof(float)*cols, NULL); 

}


void Fit::set_Npar(int par){

    this->Npar = par;

}

void Fit::set_earthquake_time(double tEq){

    this->t_Eq = tEq;

}

void Fit::set_Tau(double tau){

    this->Tau = tau;

}

void Fit::set_temporal_coherence_threshold(float tcor){

    this->tcorrThreshold = tcor;

}

void Fit::Heaviside(arma::mat &H){

     H.zeros();
     for (int ii=0; ii<nbands; ii++){
        if (acquisition_times(ii,0)>t_Eq)
            H(ii,0) = 1.0;
     }

}



void Fit::PostSeismic(arma::mat &P, arma::mat &H, double Tau){

     P.zeros();
     double pp;
     double eps = 0.000001;

     for (int ii=0; ii<nbands; ii++){
        pp = log10(1+(acquisition_times(ii,0)-t_Eq)/Tau);

        if (std::isnan(pp))
             pp = eps;

        if (std::isinf(pp))
             pp = eps;

        P(ii,0) = pp*H(ii,0);
     } 

}

void Fit::define_design_matrix(arma::mat &G, arma::mat &G_inv, str model){

  double t0 = acquisition_times(0,0);

  if (model.compare("linear") == 0){
     std::cout << "linear model" << std::endl;
     G.ones(nbands, 2);
     G.col(1) = acquisition_times.col(0) - t0;

  } else if (model.compare("linear_seasonal") == 0){
     std::cout << "linear + seasonal" << std::endl;
     G.ones(nbands, 6);
     
     G.col(1) = acquisition_times.col(0) - t0;
     G.col(2) = sin(2.0*PI*acquisition_times.col(0));
     G.col(3) = cos(2.0*PI*acquisition_times.col(0));
     G.col(4) = sin(4.0*PI*acquisition_times.col(0));
     G.col(5) = cos(4.0*PI*acquisition_times.col(0));

  } else if (model.compare("linear_coseismic") == 0){
     std::cout << "linear + coseismic" << std::endl;
     arma::mat H(nbands,1);  
     Heaviside(H);

     G.ones(nbands, 3);
     G.col(1) = acquisition_times.col(0) - t0;
     G.col(2) = H.col(0);

  } else if (model.compare("linear_coseismic_postseismic") == 0){
     std::cout << "linear + coseismic + post-seismic" << std::endl;
     std:: cout << "temporal model : a + v*(t-t_0) + H*(c + log(1+(t-tEq)/Tau)) "  << std::endl; 
     std::cout << "nbands: " << nbands << std::endl;
     arma::mat H(nbands,1);
     std::cout << "H" << std::endl;
     Heaviside(H);
     std::cout << "Heaviside" << std::endl;
     arma::mat P(nbands,1);
     std::cout << "P" << std::endl;
     PostSeismic(P, H, Tau);
     std::cout << "Post" << std::endl;
     G.ones(nbands, 4);
     G.col(1) = acquisition_times.col(0) - t0;
     G.col(2) = H.col(0);
     G.col(3) = P.col(0);
     std::cout << "G" << std::endl;
  }
  //std::cout << "rank G: " << arma::rank(G) << std::endl; 
  //G_inv = arma::pinv(G);
  bool status = arma::pinv(G_inv, G);
  std::cout << "status: " << status<< std::endl;
  std::cout << "G_inv" << std::endl;
  this->Npar = G.n_cols;
  std::cout << "Npar: " << Npar << std::endl;
}

void Fit::set_G(double *G_temp, int n_rows, int n_cols){

    //arma::mat GG(&G_temp[0], n_rows, n_cols);
    arma::mat GG(&G_temp[0], n_cols, n_rows);
    this->G = GG.t(); 
    std::cout << "G matrix " << std::endl;
    std::cout << G << std::endl;
}

void Fit::fit_timeseries(){

    arma::mat G_inv_dbl;
    arma::Mat<float> G_inv;
    arma::Mat<float> data;
    bool status = arma::pinv(G_inv_dbl, G);
    G_inv = arma::conv_to<arma::Mat<float>>::from(G_inv_dbl);

    for (int yoff=0; yoff < rows; yoff += blocksize)
    {
        std::cout << "At line " << yoff << "\n";
        int inysize = blocksize;
        if ((yoff+inysize) > rows){
            inysize = rows - yoff;
        } 
        
        get_data_block(inysize, yoff, data);
        arma::Mat<float> X = G_inv*(data.t());
        arma::Mat<float> output = arma::trans(X);
        write2datset(output, Npar, yoff, inysize); 
    
    } 

}
void Fit::fit_timeseries(str model){

    // solving for a linear system of equations data = G * X
    
    //define required matrtces
    
    //arma::mat X;
    std::cout << "define X " << std::endl;
    //arma::Mat<float> X;
    std::cout << "define G " << std::endl;
    arma::mat G;
    //arma::mat G_inv;
    std::cout << "define G_inv_dbl " << std::endl;
    arma::mat G_inv_dbl;
    std::cout << "define G_inv " << std::endl;
    arma::Mat<float> G_inv;
    std::cout << "define data " << std::endl;
    arma::Mat<float> data;
    std::cout << "define rmse " << std::endl;
    //arma::mat rmse;

    std::cout << "Design matrix " << std::endl;
    define_design_matrix(G, G_inv_dbl,  model);
    G_inv = arma::conv_to<arma::Mat<float>>::from(G_inv_dbl);

    std::cout << "Number of un-known parameters to be estimated: " << Npar << std::endl;

    for (int yoff=0; yoff < rows; yoff += blocksize)
    //for (int yoff=0; yoff < 10; yoff += blocksize)
    {
        std::cout << "At line " << yoff << "\n";
        int inysize = blocksize;
        if ((yoff+inysize) > rows){
            inysize = rows - yoff;
        }

        std::cout << "block size  " << inysize << "\n";
        
        // get a block of data
        std::cout << "getting a block of data ..." << std::endl;
        get_data_block(inysize, yoff, data);

        std::cout << "estimate un-knowns ..." << std::endl;
        std::cout << "data shape: " << data.n_rows  << " x " << data.n_cols << std::endl;
        std::cout << "G_inv shape: " << G_inv.n_rows << " x " << G_inv.n_cols << std::endl; 
        std::cout << "G shape: " << G.n_rows << " x " << G.n_cols << std::endl;
        //X = G_inv*(arma::trans(data));
        arma::Mat<float> X = G_inv*(data.t());

        std::cout << "X shape "<< X.n_rows << " x "  << X.n_cols << std:: endl;

        //rmse = arma::sqrt(arma::sum(arma::pow(arma::trans(data) - G*X, 2), 0)/nbands);
        //rmse = G*X; //arma::trans(data);
        //std::cout << " rmse: " << rmse.n_rows << " x "  << rmse.n_cols << std:: endl;

        std::cout << "writing estimated parameters ..." << std::endl;
        arma::Mat<float> output = arma::trans(X);
        write2datset(output, Npar, yoff, inysize);
        //write2datset(arma::trans(X), Npar, yoff, inysize);  
        //write2datset(arma::trans(X), rmse, Npar, yoff, inysize);
        std::cout << "writing estimated parameters finished" << std::endl;
   }
   
   std::cout << "end of fit_timeseries " << std::endl;

   /*
   std::cout << "closing the dataset" << std::endl;
   
    if (inDataset != NULL)
        GDALClose(inDataset);

    if (outDataset != NULL)
        GDALClose(outDataset);

    if (tcorDataset !=NULL)
        GDALClose(tcorDataset);

    std::cout << "Destroy the driver" << std::endl;
    GDALDestroyDriverManager();
    */
    G_inv.reset();
    G.reset();
    data.reset();
    G_inv_dbl.reset();
    

}



void Fit::fit_timeseries(str model, float &meanRMSE){

    // solving for a linear system of equations data = G * X

    // define required matrtces
    arma::mat X;
    arma::mat G;
    arma::mat G_inv;
    arma::Mat<float> data;
    arma::Mat<float> tcorr;
    arma::mat residual;
    arma::mat rmse;
    float sumRMSE=0;
    int nCoherentPixels=0;
    define_design_matrix(G, G_inv,  model);

    std::cout << "Number of un-known parameters to be estimated: " << Npar << std::endl;

    //#pragma omp parallel for\
    //    default(shared) 
    //for (int yoff=0; yoff < rows; yoff += blocksize)
        for (int yoff=3000; yoff < 3500; yoff += blocksize)
        {
            std::cout << "At line " << yoff << "\n";
            int inysize = blocksize;
            if ((yoff+inysize) > rows){
               inysize = rows - yoff;
            }

            std::cout << "block size  " << inysize << "\n";

            // get a block of data
            std::cout << "getting a block of data ..." << std::endl;
            get_data_block(inysize, yoff, data, tcorr);

            std::cout << "estimate un-knowns ..." << std::endl;
            X = G_inv*(arma::trans(data));
     
            std::cout << X.n_rows << " x "  << X.n_cols << std:: endl;
        
            // compute the residual
            rmse = arma::sqrt(arma::sum(arma::pow(arma::trans(data) - G*X, 2), 0)/nbands);
            float sum = 0;
            int nPixels = 0;

            // sum of the rmse at coherent pixels
            get_sum(arma::trans(rmse), tcorr, tcorrThreshold, sum, nPixels);
            nCoherentPixels = nCoherentPixels + nPixels;
            sumRMSE = sumRMSE + sum;

        }

    //compute average rmse at coherent pixels
    std::cout << "nCoherentPixels : " << nCoherentPixels << std::endl;
    if (nCoherentPixels>0)
        meanRMSE = sumRMSE/nCoherentPixels;
    else
        meanRMSE = 1000000.0;
    std::cout << "sumRMSE : " << std::setprecision(9) << sumRMSE << std::endl;
    std::cout << "meanRMSE : " << std::setprecision(9) << meanRMSE << std::endl;
}

void Fit::get_sum(arma::mat rmse, arma::Mat<float> tcorr, float threshold , float &sum, int &count){
     std::cout << "sum of rmse of the block" << std::endl;
     std::cout << rmse.n_rows << " x "  << rmse.n_cols << std:: endl;
     std::cout << tcorr.n_rows << " x "  << tcorr.n_cols << std:: endl;

     if (rmse.n_rows != tcorr.n_rows)
        std::cout << "the two vectors have different sizes " << std::endl;
     
     if (rmse.n_cols != tcorr.n_cols)
        std::cout << "the two vectors have different sizes " << std::endl;
     
     int nelem = rmse.n_rows;
     int num_coherent_pixels = 0;
     float sum_rmse = 0;
     for (int ii=0; ii<nelem; ii++){
        if (tcorr(ii,0) > threshold){
            sum_rmse = sum_rmse + rmse((ii,0));
            num_coherent_pixels++;
        }
     }
   
     sum = sum_rmse;
     count = num_coherent_pixels;
}


void Fit::estimate_Tau(str model){
     //str outDs = "temp_rmse";
     //bool writeFlag = false;
     double minTau = 0.03;
     double maxTau = 0.3;
     double tauStep = 0.05;
     int nr_search = (maxTau - minTau)/tauStep;

     arma::mat rmse_vs_tau(nr_search,2);
     std::cout << "searching for Tau between " << minTau << " and " << maxTau << std::endl; 
     std::cout << "Total number of search: " << nr_search << std::endl;
     double tau;
     float avg_rmse;
     for (int i=0; i<nr_search; i++){
        std::cout << "at " << i << " out of " << nr_search << std::endl;
        tau = minTau + i*tauStep; 
        Fit::set_Tau(tau);
        Fit::fit_timeseries(model, avg_rmse);
        rmse_vs_tau(i,0) = tau;
        rmse_vs_tau(i,1) = avg_rmse;
     }

     rmse_vs_tau.print();    
     arma::urowvec idx_min = arma::index_min(rmse_vs_tau); 
     std::cout << "index of minimum rmse is: " << idx_min << std::endl;
     this->Tau = minTau + idx_min(1)*tauStep;
     std::cout << "Estimated Tau : " << Tau<< std::endl;

     Fit::write_Tau(rmse_vs_tau);
}


void Fit::write_Tau(arma::mat X){

    int NumPar = X.n_cols; 
    int n_cols = 1;
    int n_ysize = X.n_rows;


    //
    GDALAllRegister();
    GDALDriver *poDriver = (GDALDriver*) GDALGetDriverByName("ENVI");
    char **mOptions = NULL;
    mOptions = CSLSetNameValue(mOptions, "INTERLEAVE", "BIP");
    mOptions = CSLSetNameValue(mOptions, "SUFFIX", "ADD");
    tauDataset = (GDALDataset*) poDriver->Create(tauDS.c_str(), n_cols, n_ysize, NumPar,
                                                GDT_Float64, mOptions);

    CSLDestroy(mOptions);
    //

    for (int ii =0; ii<NumPar; ii++){
    int status = tauDataset->GetRasterBand(ii+1)->RasterIO(GF_Write, 0, 0,
                    n_cols, n_ysize,
                    (void*) (&X.col(ii)[0]),
                    n_cols, n_ysize, GDT_Float64,
                    0 , 0, NULL);

     }
     
}



/*
     
Fit::evalauate_rmse(str datasetName, str qualityName, double &avg_rmse){

     //dataset = reinterpret_cast<GDALDataset *>( GDALOpenShared( datasetName.c_str(), GA_ReadOnly));

     //poBand->GetStatistics( false, true, &dfMin, &dfMax, nullptr, nullptr );

     //there should be an easy way to compute the average of a dataset using gdal
     // also there might be a way to do this given a  mask or a quality file with a threshold

    // for now let's compute the average over coherent pixels by looping over pixels
    
   data.set_size(cols*ysize, nbands);

   for (int bb=0; bb< nbands; bb++){
            int status = inDataset->GetRasterBand(bb+1)->RasterIO( GF_Read, 0, yoff,
                                cols, ysize,
                                (void*) (data.colptr(bb)),
                                cols, ysize, GDT_Float32,
                                sizeof(float),
                                sizeof(float)*cols, NULL);
            data.col(bb) = data.col(bb); //+ reference_phase(bb,0);
        }

}
*/
