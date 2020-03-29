// ----------------------------------------------------------------------------
// Copyright (c) 2018-, California Institute of Technology ("Caltech"). 
// U.S. Government sponsorship acknowledged.
// All rights reserved.
// 
// Author(s): Heresh Fattahi
// ----------------------------------------------------------------------------

#include "Unwrap.h"

void Unwrap::Set_connComp_dataset(str connectedComponents){

    this->connCompDS = connectedComponents;

}

void Unwrap::Set_mask_file(str maskFile){

    this->maskDS = maskFile;

}

void Unwrap::OpenDataset(){

    // Open the connected component dataset
    GDALAllRegister();
    inDataset = reinterpret_cast<GDALDataset *>( GDALOpenShared( connCompDS.c_str(), GA_ReadOnly));

    if (inDataset == NULL)
    {
        std::cout << "Cannot open stack file { " <<  connCompDS  << " } with GDAL for reading. \n";

        std::cout << "GDALOpen failed - " << connCompDS << "\n";
        std::cout << "Exiting with error code .... (102) \n";
        GDALDestroyDriverManager();
        //return 102;
    }
       
    this->cols = inDataset->GetRasterXSize();
    this->rows = inDataset->GetRasterYSize();
    this->nbands = inDataset->GetRasterCount();

    // 
    

}

void Unwrap::CreateMaskDataset(){

    str driverName = "ENVI";
    GDALDataType dtype = GDT_Byte;

    CreateMaskDataset(maskDS,  cols, rows , 1, driverName, dtype);

}

void Unwrap::CreateMaskDataset(str outputDS, int xSize, int ySize , int zSize, str driverName, GDALDataType dtype){

    GDALAllRegister();
    GDALDriver *poDriver = (GDALDriver*) GDALGetDriverByName(driverName.c_str());
    char **mOptions = NULL;
    mOptions = CSLSetNameValue(mOptions, "INTERLEAVE", "BIP");
    mOptions = CSLSetNameValue(mOptions, "SUFFIX", "ADD");
    maskDataset = (GDALDataset*) poDriver->Create(outputDS.c_str(), xSize, ySize, zSize,
                                                dtype, mOptions);
    if (maskDataset == NULL)
    {
        std::cout << "Could not create output dataset {" << outputDS << "} \n";
        std::cout << "Exiting with non-zero error code ... 104 \n";

        GDALClose(inDataset);
        GDALClose(maskDataset);
        //GDALDestroyDriverManager();
    }

    CSLDestroy(mOptions);

}

void Unwrap::CloseDataset(){
    std::cout << "Closing the datasets " << std::endl;
    if (inDataset != NULL)
        GDALClose(inDataset);

    if (maskDataset != NULL)
        GDALClose(maskDataset);

    //GDALDestroyDriverManager();

}

void Unwrap::write2datset(arma::Mat<short> X, int yoff, int ysize){

    int status = maskDataset->GetRasterBand(1)->RasterIO(GF_Write, 0, yoff,
                    cols, ysize,
                    (void*) (&X[0]),
                    cols, ysize, GDT_Byte,
                    0 , 0, NULL);

}

void Unwrap::compute_mask(int ysize, int yoff, arma::Mat<short> &mask){
    arma::Mat<short> data;
    data.set_size(cols*ysize, 1);

    for (int bb=0; bb< nbands; bb++){
        int status = inDataset->GetRasterBand(bb+1)->RasterIO( GF_Read, 0, yoff,
                                cols, ysize,
                                (void*) (data.memptr()),
                                cols, ysize, GDT_Byte,
                                sizeof(unsigned char),
                                sizeof(unsigned char)*cols, NULL);

        //mask = mask%data;
        mask.elem( find(data==0) ).zeros();
        //mask.elem( find(data==0)) = mask.elem( find(data==0)) + 1;
        }
}

void Unwrap::ComputeCommonMask(){

    arma::Mat<short> data;
    for (int yoff=0; yoff < rows; yoff += blocksize){
        std::cout << "At line " << yoff << "\n";
        int inysize = blocksize;
        if ((yoff+inysize) > rows){
            inysize = rows - yoff;
        }
    
        arma::Mat<short> mask(inysize*cols,1);
        mask.ones();
        compute_mask(inysize, yoff, mask);
        write2datset(mask, yoff, inysize); 

    }
}

