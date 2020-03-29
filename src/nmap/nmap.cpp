/**\file ampdispersion.c
 * \author Piyush Agram.
 *  */

#include "nmap.hpp"
#include "gdal.h"
#include "gdal_priv.h"
#include "cpl_conv.h"
#include <armadillo>

//Detect if compiling with nvcc 
#ifdef BUILD_WITH_CUDA
    #define RUN_CUDA_NMAP 1
#else
    #define RUN_CUDA_NMAP 0
#endif

int nmap_process(nmapOptions *opts)
{

    //Flag for test
    bool useKS = false;

    //Print user options to screen
    opts->print();

    if (opts->method.compare("KS2") == 0)
    {
        std::cout << "Using Kolmogorov-Smirnov 2-sample test\n";
        useKS = true;
    }
    else if (opts->method.compare("AD2") == 0)
    {
        std::cout << "Using Anderson-Darling 2-sample test \n";
    }
    else
    {
        std::cout << "Statistics method can be KS2 or AD2\n";
        std::cout << "Unknown method: " << opts->method << "\n";
        std::cout << "Returning with non-zero error code \n";
        return 1;
    }

    int cols, rows, nbands;
    int blockysize, boxesperblock;


    //First thing to do is to read dimenions.
    GDALDataset* inDataset = NULL;
    GDALDataset* ncountDataset = NULL;
    GDALDataset* wtsDataset = NULL;
    GDALDataset* mskDataset = NULL;

    //Workers
    KS2sample *ksworkers;
    AD2unique *adworkers;

    //Clock variables
    double t_start, t_end;

    //Register GDAL drivers
    GDALAllRegister();


    //Determine the number of uint32 bytes needed to store weights
    int nulong = ceil( ((2*opts->Ny+1)*(2*opts->Nx+1))/ 32.0);

    std::cout << "Number of uint32 bytes for mask: " << nulong << "\n";

    //Open the stack dataset
    inDataset = reinterpret_cast<GDALDataset *>( GDALOpenShared( opts->inputDS.c_str(), GA_ReadOnly));

    if (inDataset == NULL)
    {
        std::cout << "Cannot open stack file { " <<  opts->inputDS << " } with GDAL for reading. \n";

        std::cout << "GDALOpen failed - " << opts->inputDS << "\n";
        std::cout << "Exiting with error code .... (102) \n";
        GDALDestroyDriverManager();
        return 102;
    }

    cols = inDataset->GetRasterXSize();
    rows = inDataset->GetRasterYSize();
    nbands = inDataset->GetRasterCount();

    std::cout << "Number of rows  = " << rows << "\n";
    std::cout << "Number of cols  = " << cols << "\n";
    std::cout << "Number of bands = " << nbands << "\n";

    checkLongSetting();

    //Check GPU settings
    bool useGPU = !(opts->noGPU);

    if(RUN_CUDA_NMAP)
    {
        if(useGPU)
        {
            std::cout << "User has requested GPU and a device seems to be available. \n";
            std::cout << "Will execute GPU version. \n";
        }
        else
        {
            std::cout << "User had requested no GPU but device is available. \n";
            std::cout << "Executing CPU version. \n";
        }
    }
    else
    {
        std::cout << "NO GPU Available - executing CPU version \n";
        useGPU = false;
    }


    //See if mask is available
    if (!opts->maskDS.empty())
    {
        if (strcasecmp(opts->maskDS.c_str(), "None") != 0)
        {
            mskDataset = reinterpret_cast<GDALDataset *>( GDALOpenShared( opts->maskDS.c_str(), GA_ReadOnly));
            if (mskDataset == NULL)
            {
                std::cout << "Cannot open mask file { " << opts->maskDS << " } with GDAL for reading. \n";
                std::cout << "Exiting with error code .... (102) \n";
                GDALDestroyDriverManager();
                return 102;
            }

            int cc = mskDataset->GetRasterXSize();
            int rr = mskDataset->GetRasterYSize();
            int bb = mskDataset->GetRasterCount();
            int code = 0;

            if (cc != cols)
            {
                std::cout << "Mask file width does not match stack size width \n";
                code = 104;
            }

            if (rr != rows)
            {
                std::cout << "Mask file length does not match stack size length \n";
                code = 105;
            }

            if (bb != 1)
            {
                std::cout << "Mask file has more than one band \n";
                code = 106;
            }

            if (code != 0)
            {
                std::cout << "Exiting with error code .... (" << code << ")\n";
                GDALClose(inDataset);
                GDALClose(mskDataset);
                GDALDestroyDriverManager();
                return code;
            }
        }
    }


    //Determine blocksizes
    boxesperblock = int((opts->memsize * 1.0e6)/cols) / (opts->blocksize * 4 * (nbands+2+nulong)) ;
    blockysize = boxesperblock * opts->blocksize;
    if (blockysize < opts->blocksize)
    {
        blockysize = opts->blocksize;
        boxesperblock = 1;
    }

    if (blockysize > rows)
    {
        blockysize = rows;
        boxesperblock = 1;
    }

    std::cout << "Block size = " << blockysize << " lines \n";

    int totalblocks = ceil( rows / (1.0 * blockysize));
    std::cout << "Total number of blocks to process: " << totalblocks << "\n";

    //Start the clock
    t_start = getWallTime();

    //Temporary array for reading in complex data
    arma::cx_fvec cpxdata(cols*blockysize);             //To read in one SLC
    arma::Col<unsigned char> zeromask(cols*blockysize); //Mask for one SLC

    //Actual arrays
    arma::fmat amp(nbands,cols*blockysize);             //Amplitudes of all bands
    arma::Mat<int> count(cols,blockysize);              //Count for all pixels
    arma::Mat<unsigned int> wts(nulong, cols*blockysize);   //Weights for all pixels


    //Start block-by-block processing
    int blockcount = 0;
    int status;
    int nthreads = 0;

    //Store the normalizing constants
    arma::vec alpha(nbands);

    //Populate norms
    for (int bb=0; bb < nbands; bb++)
    {
        const char *bandconst = inDataset->GetRasterBand(bb+1)->GetMetadataItem( "amplitudeConstant", "slc");
        double bconst = 0.0;

        if (bandconst != NULL)
        {
             bconst = CPLScanDouble(bandconst, 20);
        }
        
        if(bconst <= 0.0)
        {
            bconst = 1.0;
            std::cout << "No calibration constant found for band " << bb+1 <<". Setting to 1.0. \n";
        }
        alpha[bb] = bconst;
    }

    double normval = alpha[0];
    for (int bb=0; bb < nbands; bb++)
    {
            alpha[bb] /= normval;
    }
    alpha[0] = 1.0;

    for(int bb=0; bb < nbands; bb++)
        std::cout << "Band " << bb+1 << ": " << alpha[bb] << "\n";

   
    //Creation of count dataset 
    GDALDriver *poDriver = (GDALDriver*) GDALGetDriverByName("ENVI");
    char **mOptions = NULL;
    mOptions = CSLSetNameValue(mOptions, "INTERLEAVE", "BIP");
    mOptions = CSLSetNameValue(mOptions, "SUFFIX", "ADD");

    ncountDataset = (GDALDataset*) poDriver->Create(opts->ncountDS.c_str(), cols, rows, 1, GDT_Int16, mOptions); 
    if (ncountDataset == NULL)
    {
        std::cout << "Could not create count dataset {" << opts->ncountDS << "} \n";
        std::cout << "Exiting with non-zero error code ... 104 \n";

        GDALClose(inDataset);
        GDALDestroyDriverManager();
        return 104;
    }

    wtsDataset = (GDALDataset*) poDriver->Create(opts->wtsDS.c_str(), cols, rows, nulong, GDT_UInt32, mOptions);
    if (wtsDataset == NULL)
    {
        std::cout << "Could not create weights dataset {" << opts->wtsDS << "}\n";
        std::cout << "Exiting with non-zero error code ... 105 \n";
        GDALClose(inDataset);
        GDALClose(ncountDataset);
        GDALDestroyDriverManager();
        return 105;
    }

    CSLDestroy(mOptions);

    nthreads = numberOfThreads();
    //Setup workers
    if (useKS)
    {
        ksworkers = new KS2sample[nthreads];
        for(int tt=0; tt < nthreads; tt++)
        {
            ksworkers[tt].length = nbands;
        }
    }
    else
    {
        adworkers = new AD2unique[nthreads];
        for(int tt=0; tt<nthreads; tt++)
        {
            adworkers[tt].init(nbands);
        }
    }

    //Useful variables
    int Ny = opts->Ny;
    int Nx = opts->Nx;
    int Wy = 2*Ny+1;
    int Wx = 2*Nx+1;
    double thresh = opts->prob;

    const Ulongmask bitmask = {Ny,Nx};

    //Block-by-block processing
    int yoff = 0;

#ifdef BUILD_NMAP_WITH_CUDA
    //Lock GPU if required
    if (useGPU) lockGPU();
#endif

    while( yoff < rows)
    {
        //Increment block counter
        blockcount++;

        //Init to zeros
        amp.zeros();
        count.zeros();
        wts.zeros();
        zeromask.ones();

        
        int inysize = blockysize;

        //Number of lines to read for the block
        if ((yoff+inysize) > rows)
            inysize = rows - yoff;

//        std::cout << "Block " << blockcount << ": " << inysize << " lines \n";

        //Read in the mask
        if (mskDataset != NULL)
        {
            status = mskDataset->GetRasterBand(1)->RasterIO( GF_Read,
                    0, yoff,
                    cols, inysize,
                    (void*) (zeromask.memptr()),
                    cols, inysize, GDT_Byte,
                    sizeof(unsigned char),
                    sizeof(unsigned char)*cols, NULL);

            if (status != 0)
            {
                std::cout << "Error reading mask band at line "<< yoff << "\n";
                std::cout << "Exiting with error code .... (108) \n";
                GDALClose(inDataset);
                GDALClose(ncountDataset);
                GDALClose(wtsDataset);
                GDALDestroyDriverManager();
                return 108;
            }
        }

        for(int bb=0; bb < nbands; bb++)
        {
            //Read in the data
            status = inDataset->GetRasterBand(bb+1)->RasterIO( GF_Read,
                        0, yoff,   
                        cols, inysize,
                        (void*) (cpxdata.memptr()), 
                        cols, inysize, GDT_CFloat32,
                        sizeof(std::complex<float>),
                        sizeof(std::complex<float>)*cols, NULL);

            if (status != 0)
            {
                std::cout << "Error reading data from band " << bb+1 << " at line " 
                    << yoff << "\n";
                std::cout << "Exiting with error code .... (108) \n";
                GDALClose(inDataset);
                GDALClose(ncountDataset);
                GDALDestroyDriverManager();
                return 108;
            }

//            std::cout << "Done reading band : " << bb+1 << "\n"; 

            //Read in normalized amplitude to amplitude array
    #pragma omp parallel for\
            default(shared)
            for(int jj=0;jj<(inysize*cols);jj++)
            {
                float val;
                val = std::abs(cpxdata(jj)) / alpha(bb);
                zeromask[jj] = (zeromask[jj] != 0) && (val != 0.) && (!std::isnan(val));
                if (zeromask[jj] != 0)
                {
                    amp(bb,jj) = val;
                }
            }
        }

        //Standard CPU-based processing
        if (!useGPU)
        {

            //Sort all the amp arrays because the Goodness-of-fit tests need this
    #pragma omp parallel for\
            default(shared)
            for(int jj=0; jj<(inysize*cols); jj++)
            {
                if(zeromask[jj] == 0) continue;

                float *ptr = amp.colptr(jj);
                std::sort(ptr, ptr+nbands);
            }
                

            //Fast computation
            //We only loop over half 2 quadrants
            //By symmetry similarity between pixel p and q
            //is same as similarity between q and p
            for(int ii=0; ii<=Ny; ii++)
            {
                int ymax = inysize-ii;
                int start = (ii==0)? 0:-Nx;
            
                for(int jj=start; jj<=Nx; jj++)
                {
                    int xmin = std::max(-jj,0);
                    int xmax = std::min(cols-jj, cols);
   
//                  std::cout << "dy = " << ii <<", dx = " << jj << "\n";
    #pragma omp parallel for\
            default(shared)
                    for(int pp=0; pp < (ymax*cols); pp++)
                    {

                        //Check zeromask
                        if (zeromask[pp] == 0) continue;

                        int cenii = pp/cols;
                        int cenjj = pp%cols;

                        if ((cenjj < xmin) || (cenjj>=xmax)) continue;

                        int refii = (cenii+ii);
                        int refjj = (cenjj+jj);
                        int qq = refii*cols + refjj;


                        //Also check zeromask here
                        if (zeromask[qq]==0) continue;

                        float* cenpix = amp.colptr(pp);
                        float* refpix = amp.colptr(qq);


                        //Always include the center pixel
                        if (pp == qq)
                        {
                            count(pp) += 1;
                            bitmask.setbit( wts.colptr(pp), ii, jj, true);
                        }
                        else
                        {
                            int threadnum = omp_get_thread_num();
                            double prob;

                            if (useKS)
                            {
                                prob = ksworkers[threadnum].test(cenpix, refpix);
                            }
                            else
                            {
                                prob = adworkers[threadnum].test(cenpix, refpix);
                            }
//                          std::cout << cenii << "   " << cenjj << "  " << prob << "   " << prob1 << "\n";

                            //Check if similar
                            if (prob >= thresh)
                            {
                                count(pp) +=1;
                                count(qq) += 1;

                                bitmask.setbit( wts.colptr(pp), ii, jj, true);
                                bitmask.setbit( wts.colptr(qq), -ii, -jj, true);
                            }
                        }
                    }
                }
            }
        }
#ifdef BUILD_NMAP_WITH_CUDA
        else    //Use GPU for processing
        {
            nmapProcessBlock(amp.memptr(), zeromask.memptr(),
                        cols, inysize, nbands,
                        count.memptr(), wts.memptr(),
                        nulong, thresh,
                        Nx, Ny);

        }
#endif
        //Block logic
        int firstlinetowrite;
        int linestowrite;
        int rollback;

        //For first block include the top half window.
        if (blockcount == 1)
        {
            firstlinetowrite = 0;
            linestowrite = inysize - Ny;
            rollback = Ny;
        
            if( (yoff+blockysize) >= rows) //There is only one block, write the whole thing
            {
                linestowrite = inysize;
                rollback = 0;
            }

        }
        else if( (yoff+blockysize) >= rows) //Last block. Write bottom half window.
        {
            firstlinetowrite = Ny;
            linestowrite = inysize-Ny;
            rollback = 0;
        }
        else
        {
            firstlinetowrite = Ny;
            linestowrite = inysize-2*Ny;
            rollback = 0;
        }


        //Print debugs
/*        std::cout << "Block number: " << blockcount << "\n";
        std::cout << "yoff: " << yoff << "\n";
        std::cout << "firstlinetowrite: " << firstlinetowrite << "\n";
        std::cout << "linestowrite: " << linestowrite << "\n";
        std::cout << "rollback: " << rollback << "\n";

        std::cout << "Gdal block \n";
        std::cout << "First line : " << yoff + firstlinetowrite << "\n";
        std::cout << "Number lines: " << linestowrite << "\n";
        std::cout << "Next read start: " << yoff + linestowrite - rollback << "\n"; */

        //Write stuff to output files


        //Write to count file
        status = ncountDataset->GetRasterBand(1)->RasterIO(GF_Write, 0, yoff+firstlinetowrite,
                    cols, linestowrite,
                    (void*) (count.colptr(firstlinetowrite)), 
                    cols, linestowrite, GDT_Int32,
                    sizeof(int), sizeof(int)*cols, NULL);

        status = wtsDataset->RasterIO(GF_Write,0,yoff+firstlinetowrite,
                            cols, linestowrite,
                            (void*)(wts.colptr(firstlinetowrite*cols)),
                            cols, linestowrite, GDT_UInt32,
                            nulong, NULL, 
                            sizeof(unsigned int)*nulong,
                            sizeof(unsigned int)*cols*nulong,
                            sizeof(unsigned int), NULL);

        if (status != 0)
        {
            std::cout << "Error writing wts data at line " << yoff << "\n";
            std::cout << "Exiting with error code .... (111) \n";
            GDALClose(inDataset);
            GDALClose(ncountDataset);
            GDALClose(wtsDataset);
            GDALDestroyDriverManager();
            return 111;
        }

//        std::cout << "Weights at (10,10): " << count(10,10) << " " 
//                  << wts(3,10*cols+10)<< "\n";
//        bitmask.print( wts.colptr(10*cols + 10));


        if ((yoff+blockysize) < rows)
        {
            yoff += linestowrite - rollback;
        }
        else
        {
            yoff=rows;
        }

        status = GDALTermProgress( yoff/(1.0*rows), NULL, NULL);
    }

#ifdef BUILD_NMAP_WITH_CUDA
    //Release GPU if required
    if (useGPU) unlockGPU();
#endif

    t_end = getWallTime();

    std::cout << "nmap processing time: " << (t_end-t_start)/60.0 << " mins \n";

    //Annotate the output datasets with window sizes for future reference
    status = ncountDataset->SetMetadataItem("HALFWINDOWX", std::to_string(Nx).c_str(), "ENVI");
    status = ncountDataset->SetMetadataItem("HALFWINDOWY", std::to_string(Ny).c_str(), "ENVI");
    status = wtsDataset->SetMetadataItem("HALFWINDOWX", std::to_string(Nx).c_str(), "ENVI");
    status = wtsDataset->SetMetadataItem("HALFWINDOWY", std::to_string(Ny).c_str(), "ENVI");

    //Clear up workers
    if (useKS)
    {
        delete [] ksworkers;
    }
    else
    {
        delete [] adworkers;
    }

    //Close the datasets
    GDALClose(inDataset);
    GDALClose(ncountDataset);
    GDALClose(wtsDataset);

    if (mskDataset != NULL)
        GDALClose(mskDataset);

    return(0);
       
};



/* Main driver*/
int main(int  argc, const char *argv[] ) {

    //Options  
    nmapOptions opts;
    int status;

    //Parse command line options
    status = opts.initFromCmdLine(argc, argv);

    if (status != 0)
    {
        std::cout << "Error parsing command line inputs \n";
        std::cout << "Exiting with non-zero return code .... \n";
        return (101);
    }


    //Execute 
    status = nmap_process(&opts);

    if (status != 0)
    {
        std::cout << "Error processing calamp \n";
        std::cout << "Exiting with non-zero return code ....\n";
        return (status);
    }

    return(0);
       
};
