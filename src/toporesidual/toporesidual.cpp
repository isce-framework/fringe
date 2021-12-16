/**\file toporesidual.cpp
 *  * \author Heresh Fattahi.
 *   *  */

#include "toporesidual.hpp"


int toporesidual_process(topoOptions *opts)
{

    //Print user options to screen
    opts->print();

    int cols, rows, nbands;
    int blockysize, boxesperblock;

    float rng, starting_range, range_spacing;
    float inc;
    float wvl;
    float Kmod, Cph; 

    //First thing to do is to read dimenions.
    GDALDataset* phaseDataset = NULL;          //Input stack of wrapped phase 
                                               //(can be wrapped phase time-series 
                                               //or wrapped interferogram network)
    GDALDataset* bperpDataset = NULL;          // Input perpendicular baseline dataset (3D)
    GDALDataset* incDataset = NULL;            // Input incidence angle map (2D)
    GDALDataset* cohDataset = NULL;            // Output temporal coherence dataset
    GDALDataset* deltazDataset = NULL;         // Output estimated DEM error

    //Clock variables
    double t_start, t_end; 
   
    //Make sure C++11 is being used
    //unsigned int - 32 bytes
    //checkLongSetting();
    
    //Register GDAL drivers
    GDALAllRegister(); 
    
    //Open the stack dataset
    phaseDataset = reinterpret_cast<GDALDataset *>( GDALOpenShared( opts->phaseDS.c_str(), GA_ReadOnly));
    if (phaseDataset == NULL)
    {
        std::cout << "Cannot open stack file { " <<  opts->phaseDS << " } with GDAL for reading. \n";

        std::cout << "GDALOpen failed - " << opts->phaseDS << "\n";
        std::cout << "Exiting with error code .... (102) \n";
        GDALDestroyDriverManager();
        return 102;
    }

    cols = phaseDataset->GetRasterXSize();
    rows = phaseDataset->GetRasterYSize();
    nbands = phaseDataset->GetRasterCount();

    std::cout << "Number of rows  = " << rows << "\n";
    std::cout << "Number of cols  = " << cols << "\n";
    std::cout << "Number of bands = " << nbands << "\n";    

    // Open perpendicular baseline dataset
    bperpDataset = reinterpret_cast<GDALDataset *>( GDALOpenShared( opts->bperpDS.c_str(), GA_ReadOnly));
    if (bperpDataset == NULL)
    {
        std::cout << "Cannot open stack file { " <<  opts->bperpDS << " } with GDAL for reading. \n";

        std::cout << "GDALOpen failed - " << opts->bperpDS << "\n";
        std::cout << "Exiting with error code .... (102) \n";
        GDALDestroyDriverManager();
        return 103;
    } 

    // Open incidence angle dataset    
    incDataset = reinterpret_cast<GDALDataset *>( GDALOpenShared( opts->incDS.c_str(), GA_ReadOnly));
    if (bperpDataset == NULL)
    {
        std::cout << "Cannot open stack file { " <<  opts->incDS << " } with GDAL for reading. \n";

        std::cout << "GDALOpen failed - " << opts->incDS << "\n";
        std::cout << "Exiting with error code .... (102) \n";
        GDALDestroyDriverManager();
        return 104;
    }

    // check for consistency between among input datasets


    //Determine blocksizes
    //Start the clock
    t_start = getWallTime();    
    
    // create output datasets
    {
            GDALDriver *poDriver = (GDALDriver*) GDALGetDriverByName("ENVI");
            char **mOptions = NULL;
            mOptions = CSLSetNameValue(mOptions, "INTERLEAVE", "BIP");
            mOptions = CSLSetNameValue(mOptions, "SUFFIX", "ADD");

            std::string fname = CPLFormFilename(opts->outputFolder.c_str(), opts->deltazDS.c_str(), NULL);

            std::cout << "Corr: " << fname <<  " " << fname.size() << "\n";
            deltazDataset = (GDALDataset*) poDriver->Create(fname.c_str(), cols, rows, 1, GDT_Float32, mOptions);

            if (deltazDataset == NULL)
            {
                std::cout << "Could not create temporal correlation file: " << fname << "\n";
                std::cout << "Exiting with non-zero error code ... 117 \n";

                GDALClose(phaseDataset);
                GDALClose(bperpDataset);

                GDALDestroyDriverManager();
                return 117;
            }

            std::cout << "Created : " << fname << " " << deltazDataset << "\n";

            // create output temporal coherence
            fname = CPLFormFilename(opts->outputFolder.c_str(), opts->cohDS.c_str(), NULL);
            std::cout << "Corr: " << fname <<  " " << fname.size() << "\n";
            cohDataset = (GDALDataset*) poDriver->Create(fname.c_str(), cols, rows, 1, GDT_Float32, mOptions);
            
            if (deltazDataset == NULL)
            {
                std::cout << "Could not create temporal coherence file: " << fname << "\n";
                std::cout << "Exiting with non-zero error code ... 117 \n";

                GDALClose(phaseDataset);
                GDALClose(bperpDataset);

                GDALDestroyDriverManager();
                return 117;
            }

            std::cout << "Created : " << fname << " " << cohDataset << "\n";

            CSLDestroy(mOptions);
        }    


    topofit worker(20.0, 0.1, wvl, nbands);

    //read data block by block
    int yoff = 0;
    int blockcount = 0;
    blockysize = int((opts->memsize * 1.0e6)/(cols * 8 * 6) );
    std::cout << "Computed block size based on memory size = " << blockysize << " lines \n";
    int totalblocks = ceil( rows / (1.0 * blockysize));
    std::cout << "Total number of blocks to process: " << totalblocks << "\n";

    // Array for reading complex interferogram (or time-series) stack
    arma::cx_fmat cpxdata(cols*blockysize, nbands); 
    arma::Mat<float> bperp(cols*blockysize, nbands);
    arma::Mat<float> incidenceAngle(cols*blockysize, 1);
    arma::Mat<float> coherence(cols*blockysize, 1);
    arma::Mat<float> dz(cols*blockysize, 1);
    
    int status;

    while( yoff < rows)
    {
        //Increment block counter
        blockcount++;   
        int inysize = blockysize;

        //Number of lines to read for the block
        if ((yoff+inysize) > rows)
            inysize = rows - yoff;

        std::cout << "inysize: " << inysize << std::endl;
        
        //read the block of complex data and bperp
        for(int bb=0; bb < nbands; bb++)
        {
            status = phaseDataset->GetRasterBand(bb+1)->RasterIO( GF_Read,
                        0, yoff,
                        cols, inysize,
                        (void*) (cpxdata.colptr(bb)),
                        cols, inysize, GDT_CFloat32,
                        sizeof(std::complex<float>),
                        sizeof(std::complex<float>)*cols, NULL);
             if (status != 0)
            {
                std::cout << "Error reading data from band " << bb+1 << " at line "
                    << yoff << "\n";
                std::cout << "Exiting with error code .... (108) \n";
                GDALClose(phaseDataset);
                GDALClose(bperpDataset);
                GDALClose(cohDataset);
                GDALDestroyDriverManager();
                return 108;
            }

            status = bperpDataset->GetRasterBand(bb+1)->RasterIO( GF_Read,
                        0, yoff,
                        cols, inysize,
                        (void*) (bperp.colptr(bb)),
                        cols, inysize, GDT_Float32,
                        sizeof(float),
                        sizeof(float)*cols, NULL);

        }
        status = incDataset->GetRasterBand(1)->RasterIO( GF_Read,
                        0, yoff,
                        cols, inysize,
                        (void*) (incidenceAngle.memptr()),
                        cols, inysize, GDT_Float32,
                        sizeof(float),
                        sizeof(float)*cols, NULL);

        float rng = 699993.605;
        float wvl = 0.031228381041666666;
        float Kmod, Cph;
        std::cout << "starting range: " << rng << std::endl;
        std::cout << "wavelength: " << wvl << std::endl;
        arma::cx_fvec cph(nbands);
        arma::cx_fvec resid(nbands);
        arma::fvec bb(nbands,1); 

        /*cph = cpxdata.row(100).t();
        bb = bperp.row(100).t();
        
        topofit worker1(20.0, 0.1, wvl, nbands);        
        float cc  = worker1.fit(cph.memptr(), bb.memptr(), rng, incidenceAngle(100), resid.memptr(), dzz, Cph);

        std::cout << " dzz: " << dzz
                << " Cph: " << Cph << " coh: " << cc << std::endl;
        */
        float cc; 
        topofit worker(20.0, 0.1, wvl, nbands);
        //#pragma omp parallel for\
            //default(shared)
            for(int jj=0; jj<(inysize*cols); jj++)
            {
              cph = cpxdata.row(jj).t();
              bb = bperp.row(jj).t();
              cc  = worker.fit(cph.memptr(), bb.memptr(), rng, incidenceAngle(jj), resid.memptr(), Kmod, Cph);
              dz(jj) = Kmod*(wvl*rng*sin(incidenceAngle(jj)*M_PI/180.0)/4.0/M_PI);
              coherence(jj) = cc;
            }


        status = cohDataset->GetRasterBand(1)->RasterIO(GF_Write, 0, yoff,
                    cols, inysize,
                    (void*) (&coherence[0]),
                    cols, inysize, GDT_Float32,
                    sizeof(float), sizeof(float)*cols, NULL);

        status = deltazDataset->GetRasterBand(1)->RasterIO(GF_Write, 0, yoff,
                    cols, inysize,
                    (void*) (&dz[0]),
                    cols, inysize, GDT_Float32,
                    sizeof(float), sizeof(float)*cols, NULL);

        if ((yoff+inysize) < rows)
        {
            yoff += inysize;
        }
        else
        {
            yoff=rows;
        }

     }

     return(0);

}

int main(int  argc, const char *argv[] ) {
    //Options
    topoOptions opts;
    int status;

    //Parse command line options
    status = opts.initFromCmdLine(argc, argv);
    
    if (status != 0)  //Help message
    {
        std::cout << "Error processing command line for toporesidual \n";
        std::cout << "Exiting with non-zero return code .... \n";
        return(status);
    }

    //Execute
    status = toporesidual_process(&opts);

    if (status != 0)
    {
        std::cout << "Error processing toporesidual \n";
        std::cout << "Exiting with non-zero return code ....\n";
        return (status);
    }

    return(0);

}; 


